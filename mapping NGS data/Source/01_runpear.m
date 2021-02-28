filesPath = '../Data/Samples/';
outPath = '../Data/Assembled/';


if ~isdir(outPath)
	mkdir(outPath);
end

shift_primers={''};
f_homology=''; %%% insert R primer here
r_homology=''; %%% insert F primer homology (reverse complement)
adapter_str='';
for i=1:length(shift_primers)
	adapter_str=[adapter_str sprintf('-a %s%s...%s ',shift_primers{i},f_homology,r_homology)];
end


samplesFiles = dir([filesPath '*.fastq*']);
samplesFiles = {samplesFiles.name};


samplesNames = cellfun(@(x) x(1:regexp(x,'_S[0-9]{1,3}_L00[1-4]_R1_001')-1),samplesFiles,'un',0);
samplesNames = unique(samplesNames);


for i = 1:numel(samplesNames)
	files = dir([filesPath samplesNames{i} '_*']);
	files = {files.name};
	if isempty(files)
		continue
	end
	lanes = cellfun(@(x) x(1:regexp(x,'_R[0-9]_001')-1),files,'un',0);
	lanes=unique(lanes);
    currPath = outPath;
	for f = 1:numel(lanes)
		currlane = lanes{f};
		outFile = [currPath currlane];
		inFile = [filesPath currlane];
		
		str = sprintf('bsub -q new-short -u /dev/null ''pear -f "%s" -r "%s" -o "%s"''',[inFile '_R1_001.fastq'],[inFile '_R2_001.fastq'],outFile);
		system(str);
	end
end