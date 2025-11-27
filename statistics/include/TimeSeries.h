#ifndef TIME_SERIES_H
#define TIME_SERIES_H

#include <vector>

namespace Stats {

class TimeSeries {
public:
    // Autocorrelation Function (ACF)
    // Returns autocorrelation at lag k
    static double autocorrelation(const std::vector<double>& data, int lag);
    
    // Calculate ACF for lags 0 to maxLag
    static std::vector<double> calculateACF(const std::vector<double>& data, int maxLag);
};

} // namespace Stats

#endif // TIME_SERIES_H
