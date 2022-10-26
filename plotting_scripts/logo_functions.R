convert_pwm_to_ppm <- function(pwm){
    pwm = pwm[parameter %like% 'motif']
    # convert coeffs to prob
    pwm[, e_coefficient := exp(coefficient)]
    pwm[, prob := e_coefficient/(1+e_coefficient)]
    pwm[, ppm := prob/sum(prob), by = parameter]

    return(pwm)
}

convert_ppm_to_matrix <- function(ppm){
    positions = get_positions()
    mat = matrix(0, nrow = 4, ncol = length(positions)) 
    rownames(mat) = c('A', 'T', 'C', 'G')
    colnames(mat) = positions

    for (b in rownames(mat)){
        for (pos in colnames(mat)){
            entry = ppm[base == b & parameter == pos]$ppm
            mat[b,pos] = entry
        }
    }
    return(mat)
}
