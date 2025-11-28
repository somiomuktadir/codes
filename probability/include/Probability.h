#ifndef PROBABILITY_H
#define PROBABILITY_H

#include <vector>

namespace Probability {
    // Basic Probability
    double probability(double favorable_outcomes, double total_outcomes);
    double complementaryProbability(double prob);

    // Conditional Probability
    double conditionalProbability(double prob_intersection, double prob_given);

    // Bayes' Theorem
    // P(A|B) = (P(B|A) * P(A)) / P(B)
    double bayesTheorem(double prob_B_given_A, double prob_A, double prob_B);

    // Set Operations (for events)
    // P(A U B) = P(A) + P(B) - P(A n B)
    double unionProbability(double prob_A, double prob_B, double prob_intersection);
    
    // Independent Events
    // P(A n B) = P(A) * P(B)
    double intersectionIndependent(double pA, double pB);
    double intersectionDependent(double pA, double pBgivenA); // New feature
    double totalProbability(const std::vector<double>& priors, const std::vector<double>& conditionals); // New feature

    // Random Variables
    double expectedValue(const std::vector<double>& values, const std::vector<double>& probabilities);
    double variance(const std::vector<double>& values, const std::vector<double>& probabilities);

    // Explanations
    void printFormulas();
}

#endif // PROBABILITY_H
