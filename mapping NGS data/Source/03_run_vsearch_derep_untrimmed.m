filesPath = '../Data/Trimmed_indexes_3/';
outPath = '../Data/Splicing_lib/Out_derep_untrimmed/';


if ~isdir(outPath)
	mkdir(outPath);
end

samplesDirs = dir([filesPath]);
samplesDirs(1:2)=[];
samplesDirs = {samplesDirs.name};

for i = 1:numel(samplesDirs)

	currSamplePath = [filesPath samplesDirs{i} '/'];
	currOutPath=[outPath samplesDirs{i} '/'];
	if ~isdir(currOutPath)
		mkdir(currOutPath);
	end
	
	files=dir([currSamplePath '*.fastq']);
	files={files.name};
	for j=1:length(files)
		outFile=files{j}(1:regexp(files{j},'.trim.fastq')-1);
	
		cmd = sprintf('bsub -q molgen -R "rusage[mem=4096]" ''vsearch --derep_prefix %s%s --output %s%s.fasta -minseqlen 12 -maxseqlen 300 -sizeout''',currSamplePath,files{j},currOutPath,outFile);
		system(cmd);
	end
end