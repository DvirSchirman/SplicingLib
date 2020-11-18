samplePath = '../Data/Out_derep_untrimmed/';
outPath='../Data/Out_mapped_seq/';


if ~isdir(outPath)
	mkdir(outPath);
end


samplesDirs = dir([samplePath]);
samplesDirs(1:2)=[];
samplesDirs = {samplesDirs.name};
% samplesDirs = {'ND-anc'};



for i=1:length(samplesDirs)
% for i=4:4

	currSamplePath = [samplePath samplesDirs{i} '/'];
	currOutPath=[outPath samplesDirs{i} '/'];
	if ~isdir(currOutPath)
		mkdir(currOutPath);
	end
	
	fileNames = dir([currSamplePath '*.fasta']);
	fileNames={fileNames.name};


	for j=1:length(fileNames)
		cmd = sprintf('bsub -q new-medium -R "rusage[mem=16384]" matlab -nojvm -nodesktop -nosplash -r ''map_uniques_2db_seq_SP %s %s %s''',currSamplePath, fileNames{j}, currOutPath);
		system(cmd);
	end
end
