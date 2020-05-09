function process_transcriptome(input_count_matrix, outdir)

    diary([tempname '.txt']);
    diary on;
    fprintf('Writing log to temporary location: %s\n', get(0,'DiaryFile'));

    tic;

    transcriptome = Transcriptome(input_count_matrix, struct, outdir);
    transcriptome.prepare_CARLIN_input(outdir);
    save(sprintf('%s/Transcriptome.mat', outdir), 'transcriptome');
    
    fprintf('Pipeline completed in %g seconds\n', toc);
    
    diary off;
    copyfile(get(0,'DiaryFile'), [outdir '/Log.txt']);

end