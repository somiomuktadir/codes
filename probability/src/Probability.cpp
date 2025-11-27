#include "Probability.h"
#include <iostream>
#include <stdexcept>

namespace Probability {

    double probability(double favorable_outcomes, double total_outcomes) {
        if (total_outcomes <= 0) throw std::invalid_argument("Total outcomes must be positive.");
        if (favorable_outcomes < 0) throw std::invalid_argument("Favorable outcomes cannot be negative.");
        return favorable_outcomes / total_outcomes;
    }

    double complementaryProbability(double prob) {
        if (prob < 0 || prob > 1) throw std::invalid_argument("Probability must be between 0 and 1.");
        return 1.0 - prob;
    }

    double conditionalProbability(double prob_intersection, double prob_given) {
        if (prob_given <= 0 || prob_given > 1) throw std::invalid_argument("P(Given) must be in (0, 1].");
        if (prob_intersection < 0 || prob_intersection > 1) throw std::invalid_argument("P(Intersection) must be in [0, 1].");
        if (prob_intersection > prob_given) throw std::invalid_argument("P(Intersection) cannot be greater than P(Given).");
        return prob_intersection / prob_given;
    }

    double bayesTheorem(double prob_B_given_A, double prob_A, double prob_B) {
        if (prob_B <= 0) throw std::invalid_argument("P(B) must be positive.");
        // P(A|B) = P(B|A) * P(A) / P(B)
        return (prob_B_given_A * prob_A) / prob_B;
    }

    double unionProbability(double prob_A, double prob_B, double prob_intersection) {
        if (prob_A < 0 || prob_A > 1 || prob_B < 0 || prob_B > 1 || prob_intersection < 0 || prob_intersection > 1) {
            throw std::invalid_argument("Probabilities must be between 0 and 1.");
        }
        double result = prob_A + prob_B - prob_intersection;
        if (result > 1.0 + 1e-9) throw std::invalid_argument("Resulting probability > 1. Check input values.");
        return (result > 1.0) ? 1.0 : result;
    }

    double intersectionIndependent(double prob_A, double prob_B) {
        if (prob_A < 0 || prob_A > 1 || prob_B < 0 || prob_B > 1) {
            throw std::invalid_argument("Probabilities must be between 0 and 1.");
        }
        return prob_A * prob_B;
    }

    double expectedValue(const std::vector<double>& values, const std::vector<double>& probabilities) {
        if (values.size() != probabilities.size()) throw std::invalid_argument("Values and probabilities must have same size.");
        double ev = 0.0;
        double sum_p = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            if (probabilities[i] < 0 || probabilities[i] > 1) throw std::invalid_argument("Probabilities must be in [0, 1].");
            ev += values[i] * probabilities[i];
            sum_p += probabilities[i];
        }
        if (std::abs(sum_p - 1.0) > 1e-9) throw std::invalid_argument("Probabilities must sum to 1.");
        return ev;
    }

    double variance(const std::vector<double>& values, const std::vector<double>& probabilities) {
        double ev = expectedValue(values, probabilities);
        double ev_sq = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            ev_sq += (values[i] * values[i]) * probabilities[i];
        }
        return ev_sq - (ev * ev);
    }

    void printFormulas() {
        std::cout << "\n--- Probability Formulas ---\n";
        std::cout << "Basic Probability: P(E) = Favorable / Total\n";
        std::cout << "Complementary: P(E') = 1 - P(E)\n";
        std::cout << "Conditional: P(A|B) = P(A n B) / P(B)\n";
        std::cout << "Bayes' Theorem: P(A|B) = (P(B|A) * P(A)) / P(B)\n";
        std::cout << "Union (General): P(A U B) = P(A) + P(B) - P(A n B)\n";
        std::cout << "Intersection (Independent): P(A n B) = P(A) * P(B)\n";
        std::cout << "Expected Value E[X]: sum(x_i * p_i)\n";
        std::cout << "Variance Var(X): E[X^2] - (E[X])^2\n";
        std::cout << "----------------------------\n";
    }
}
