#ifndef TEMPERATURE_PROFILE
    #define TEMPERATURE_PROFILE 1
#endif

double* TempVector;
double* TimeVector;
int nPoints;

/**
 * @brief Reads temperature profile data from given arrays.
 *
 * This function reads temperature and time data from the provided arrays
 * and stores them in dynamically allocated vectors.
 *
 * @param x Pointer to an array of time data.
 * @param y Pointer to an array of temperature data.
 * @return void
 * @note The arrays x and y must have the same length.
 */
void TemperatureProfile_Set(const double*x, const double* y) {

    if (sizeof(x) != sizeof(y)) {
        printf("Error: The arrays Time and Temperature must have the same length.\n");
        return;
    }

    nPoints = sizeof(x)/sizeof(double);
    TempVector = (double*)malloc(nPoints*sizeof(double));
    TimeVector = (double*)malloc(nPoints*sizeof(double));
    for (int i=0; i<nPoints; i++) {
        TempVector[i] = y[i];
        TimeVector[i] = x[i];
    }
}

/**
 * @brief Retrieves the temperature at a given time using linear interpolation.
 *
 * This function takes a time value and returns the corresponding temperature
 * by performing linear interpolation between the known temperature points.
 * If the time is out of the bounds of the known time points, it returns an error
 * or the last known temperature value.
 *
 * @param time The time at which the temperature is to be retrieved.
 * @return The interpolated temperature at the given time. If the time is out of bounds,
 *         it returns -1 and prints an error message.
 */
double TemperatureProfile_GetT(double time) {

    if (time < TimeVector[0]) {
        printf("Error: Time is out of bound.\n");
        return -1;
    }

    // If time is greater than the last time point, return the last temperature value.
    if (time > TimeVector[nPoints-1]) {
        return TempVector[nPoints-1];
    }

    for (int i=0; i<nPoints-1; i++) {
        if (time >= TimeVector[i] && time <= TimeVector[i+1]) {
            return TempVector[i] + (TempVector[i+1]-TempVector[i])/(TimeVector[i+1]-TimeVector[i])*(time-TimeVector[i]);
        }
    }
    return -1;
}