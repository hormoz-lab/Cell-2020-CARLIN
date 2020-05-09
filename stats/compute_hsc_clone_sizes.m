function out = compute_hsc_clone_sizes(dat)
    
    out.test_alleles = union(dat.sig_alleles.hsc_only_allele, dat.sig_alleles.hsc_derived_allele);    
    out.N_transcripts = sum(dat.summary.allele_freqs(out.test_alleles));
    out.clone_size_freqs = accumarray(dat.summary.allele_freqs(out.test_alleles), 1);
    out.transcripts_per_clone_size = out.clone_size_freqs.*[1:length(out.clone_size_freqs)]';
    assert(sum(out.transcripts_per_clone_size) ==out.N_transcripts);    
    
end