#TODO add 5' or 3' designation...for orientation

get_motif_context <- function(whole_gene_seq, trimmed_gene_seq, trim_length){
    whole_gene_seq = DNAString(whole_gene_seq)
    trimmed_gene_seq = DNAString(trimmed_gene_seq)
    trimmed_length = nchar(trimmed_gene_seq)
    original_trimmed_length = trimmed_length

    if (nchar(trimmed_gene_seq) < LEFT_NUC_MOTIF_COUNT){
        trimmed_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq)-trim_length)
        trimmed_length = nchar(trimmed_gene_seq)
    }

    left_nuc_motif = substring(trimmed_gene_seq, trimmed_length - (LEFT_NUC_MOTIF_COUNT - 1), trimmed_length) 

    possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, trim_length)
    whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
    seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-2*trim_length + 1, nchar(whole_gene_and_pnucs))

    if (nchar(seq_right_of_trim) < RIGHT_NUC_MOTIF_COUNT){
        missing_nucs = DNAString(strrep('-', RIGHT_NUC_MOTIF_COUNT - nchar(seq_right_of_trim)))
        seq_right_of_trim = c(unlist(seq_right_of_trim), unlist(missing_nucs))
    }

    right_nuc_motif = substring(seq_right_of_trim, 1, RIGHT_NUC_MOTIF_COUNT)

    motif = c(unlist(left_nuc_motif), unlist(right_nuc_motif))
    return(motif)
}


