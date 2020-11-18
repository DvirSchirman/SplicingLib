filesPath = '../Data/Assembled/';
trimPath = '../Data/Trimmed_indexes_3/';
concatPath = '../Data/samples_concat/';
f_homology='AATACGAGGCACTTACTCCG'; %%% splicing

if ~isdir(trimPath)
	mkdir(trimPath);
end
if ~isdir(concatPath)
	mkdir(concatPath);
end

indexes={'AACATCTA';'TGTTGGGA';'AAGCCATG';'GCTAAAGA'}
adapter_3='TCAATGTGACTGCGTTCCAC';


samplesFiles = dir([filesPath '*.assembled.fastq*']);
disp([filesPath '*.assembled.fastq*']);
samplesFiles = {samplesFiles.name};

samplesNames = cellfun(@(x) x(1:regexp(x,'_S[0-9]{1,3}_L00[1-4]')-1),samplesFiles,'un',0);
samplesNames = unique(samplesNames);
disp(samplesNames)


for i = 1:numel(samplesNames)
	samplesFiles=[filesPath samplesNames{i} '_*.assembled*'];
	cmd=sprintf('cat %s > %s%s.fastq',samplesFiles,concatPath,samplesNames{i});
	system(char(cmd));
	
	if ~isdir([trimPath samplesNames{i}])
        mkdir([trimPath samplesNames{i}]);
    end
    
	currPath = [trimPath samplesNames{i} '/'];
		
	inFile = sprintf('%s%s.fastq',concatPath,samplesNames{i});
		
	for j=1:length(indexes)
		outFile = sprintf('%s%s.idx%d.trim.fastq',currPath,samplesNames{i},j);
		str = sprintf('bsub -q pilpel -u /dev/null -R "rusage[mem=2048]" ''cutadapt -g %s%s -a %s$ --discard-untrimmed -o %s %s''', indexes{j},f_homology, adapter_3, outFile, inFile);
		system(str);
	end

		
end