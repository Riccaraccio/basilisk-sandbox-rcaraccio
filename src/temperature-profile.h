#ifndef TEMPERATURE_PROFILE
    #define TEMPERATURE_PROFILE 1
#endif

double* TempVector = NULL;
double* TimeVector = NULL;
int nPoints;
extern double TG0;

/**
 * @brief Reads temperature profile data from given arrays.
 *
 * This function reads temperature and time data from the provided arrays
 * and stores them in dynamically allocated vectors.
 *
 * @param x Pointer to an array of time data.
 * @param y Pointer to an array of temperature data.
 * @param size The number of data points in the arrays.
 * @return void
 * @note The arrays x and y must have the same length.
 */
void TemperatureProfile_Set(const double*x, const double* y, const int size) {

    if (y[0] != TG0 || TG0 == 0) {
        fprintf(stderr,"Error: The first temperature value must be equal to TG0.\n");
        return;
    }

    if (sizeof(x) != sizeof(y)) {
        fprintf(stderr,"Error: The arrays Time and Temperature must have the same length.\n");
        return;
    }

    nPoints = size; 
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
 * @return The interpolated temperature at the given time. If the time greater than
 * the last time, returns the last temperature.
 */
double TemperatureProfile_GetT(double time) {

    if (time < TimeVector[0]) {
        fprintf(stderr,"Error: Time is out of bound.\n");
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

/**
 * @brief Frees the memory allocated for temperature and time vectors.
 *
 * This function releases the dynamically allocated memory for the 
 * temperature vector (TempVector) and the time vector (TimeVector).
 * It should be called to avoid memory leaks when these vectors are 
 * no longer needed.
 */
void TemperatureProfile_Free() {
    free(TempVector), TempVector = NULL;
    free(TimeVector), TimeVector = NULL;
}

/**
 * @brief Prints the temperature profile data to the standard error stream.
 *
 * This function prints the temperature profile data to the standard error stream.
 * It is useful for debugging purposes to verify that the data has been read correctly.
 */
void TemperatureProfile_Print() {
    for (int i=0; i<nPoints; i++) {
        fprintf(stderr,"%g %g\n", TimeVector[i], TempVector[i]);
    }
}