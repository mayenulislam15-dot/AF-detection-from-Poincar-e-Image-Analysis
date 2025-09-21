function entropyValue = calculateEntropy(matrix)
    % Normalize to probabilities
    probabilities = matrix / sum(matrix(:));

    % Exclude zero probabilities to avoid NaN in the logarithm
    probabilities = probabilities(probabilities > 0);

    % Calculate entropy
    entropyValue = -sum(probabilities .* log2(probabilities));
end