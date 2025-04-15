//
// Created by Tsail on 6/16/2023.
//

#include "SailorMath.h"

double BinomialCoefficient(double a, int k){
    ///Returns the evaluation of the given binomial coefficient (a k) = (a*(a-1)*...*(a-k+1)/k!
    double value = 1.0;
    for (int i=0; i<k; i++){
        value *= (a-i)/(k-i);
    }
    return value;
}