function p = kde_likelihood(posterior_state, prior_state, kde_pdf_function)

adjusted_state = posterior_state - prior_state;  % Center measurement around particle state
p = kde_pdf_function(adjusted_state);  % Evaluate KDE at adjusted state

end