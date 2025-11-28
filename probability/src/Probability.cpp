#include "Probability.h"
#include "Utils.h"
#include <iostream>
#include <stdexcept>

namespace Probability {

    // Helper function to validate probability values
    bool isValidProbability(double p) {
        return p >= 0.0 && p <= 1.0;
    }

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
        Utils::log("Calculating Bayes' Theorem");
        Utils::logStep("Formula: P(A|B) = P(B|A) * P(A) / P(B)");
        double result = (prob_B_given_A * prob_A) / prob_B;
        Utils::logStep("Result: " + std::to_string(result));
        return result;
    }

    double unionProbability(double pA, double pB, double pIntersection) {
        if (!isValidProbability(pA) || !isValidProbability(pB) || !isValidProbability(pIntersection))
            throw std::invalid_argument("Invalid probability values");
        Utils::log("Calculating Union Probability P(A U B)");
        Utils::logStep("Formula: P(A) + P(B) - P(A n B)");
        Utils::logStep("Substitution: " + std::to_string(pA) + " + " + std::to_string(pB) + " - " + std::to_string(pIntersection));
        double result = pA + pB - pIntersection;
        Utils::logStep("Result: " + std::to_string(result));
        // The original code had a check for result > 1.0 + 1e-9 and clamped it.
        // The provided snippet removes this, so we follow the snippet.
        return result;
    }

    double intersectionIndependent(double pA, double pB) {
        if (!isValidProbability(pA) || !isValidProbability(pB))
            throw std::invalid_argument("Invalid probability values");
        Utils::log("Calculating Intersection (Independent) P(A n B)");
        Utils::logStep("Formula: P(A) * P(B)");
        Utils::logStep("Substitution: " + std::to_string(pA) + " * " + std::to_string(pB));
        double result = pA * pB;
        Utils::logStep("Result: " + std::to_string(result));
        return result;
    }

    double intersectionDependent(double pA, double pBgivenA) {
        if (!isValidProbability(pA) || !isValidProbability(pBgivenA))
            throw std::invalid_argument("Invalid probability values");
        Utils::log("Calculating Intersection (Dependent) P(A n B)");
        Utils::logStep("Formula: P(A) * P(B|A)");
        Utils::logStep("Substitution: " + std::to_string(pA) + " * " + std::to_string(pBgivenA));
        double result = pA * pBgivenA;
        Utils::logStep("Result: " + std::to_string(result));
        return result;
    }

    double totalProbability(const std::vector<double>& priors, const std::vector<double>& conditionals) {
        if (priors.size() != conditionals.size())
            throw std::invalid_argument("Size mismatch between priors and conditionals");
        
        Utils::log("Calculating Total Probability");
        Utils::logStep("Formula: Sum(P(Ai) * P(B|Ai))");
        
        double total = 0.0;
        for (size_t i = 0; i < priors.size(); ++i) {
            if (!isValidProbability(priors[i]) || !isValidProbability(conditionals[i]))
                throw std::invalid_argument("Invalid probability values");
            
            double term = priors[i] * conditionals[i];
            if (Utils::isVerbose()) {
                Utils::logStep("Term " + std::to_string(i+1) + ": " + std::to_string(priors[i]) + " * " + std::to_string(conditionals[i]) + " = " + std::to_string(term));
            }
            total += term;
        }
        Utils::logStep("Total Result: " + std::to_string(total));
        return total;
    }

    double expectedValue(const std::vector<double>& values, const std::vector<double>& probabilities) {
        Utils::log("Calculating Expected Value E[X]");
        Utils::logStep("Formula: Sum(x_i * p_i)");
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
