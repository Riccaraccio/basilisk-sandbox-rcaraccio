#include "temperature-profile.h"

int main () {
    double x[] = {0, 1, 2, 3, 4};
    double y[] = {20, 25, 30, 35, 40};
    TemperatureProfile_Set(x, y);
    printf("Temperature at time 2.5: %f\n", TemperatureProfile_GetT(2.5));
    printf("Temperature at time 5: %f\n", TemperatureProfile_GetT(5));
    return 0;
}