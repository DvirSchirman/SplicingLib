samplePath = '../Data/Out_mapped_seq/';
outPath='../Data/Out_double_sp_eff_junctions/';

if ~isdir(outPath)
	mkdir(outPath);
end

samplesDirs = dir([samplePath]);
samplesDirs(1:2)=[];
samplesDirs = {samplesDirs.name};

for i=1:length(samplesDirs)
	currSamplePath = [samplePath samplesDirs{i} '/'];
	currOutPath=[outPath samplesDirs{i} '/'];
	if ~isdir(currOutPath)
		mkdir(currOutPath);
	end
	
	fileNames = dir([currSamplePath '*.mat']);
	fileNames={fileNames.name};

	for j=1:length(fileNames)
		cmd = sprintf('bsub -q new-short -R "rusage[mem=4096]" matlab -nojvm -nodesktop -nosplash -r ''get_double_sp_eff_junctions %s %s %s''',currSamplePath, fileNames{j}, currOutPath);
		system(cmd);
	end
end
