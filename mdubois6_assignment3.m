RBI_5UTR = 'GCUCAGUUGCCGGGCGGGGGAGGGCGCGUCCGGUUUUUCUCAGGGGACGUUGAAAUUAUUUUUGUAACGGGAGUCGGGAGAGGACGGGGCGUGCCCCGACGUGCGCGCGCGUCGUCCUCCCCGGCGCUC';

wild_type_1 = '(((...((((..((((((((..(((((((((.............))))).............(((...(((.((((((.((((((((.((((((..........)))))).))))))))))))))..)))...)))...))))..))))))))))))...)))...';
wild_type_2 =  '(((....(.((((((((((((((((((((((.............))))))......................((((((.((((((((.((((((..........)))))).)))))))))))))))).)))))..((((....)))))))))).))))..)))...';
wild_type_3 = '........(((.(((((((((((.(((((((.............)))))).)....................((((((.((((((((.((((((..........)))))).))))))))))))))...)))))....)))))))))...((((.......))))..';

RBI_5UTR_dot_bracket = rnafold(RBI_5UTR, 'MinLoopSize',7); %Step 1: Illustrate the secondary structure of RBI's 5'UTR
RBI_5UTR_diagram = rnaplot(RBI_5UTR_dot_bracket, 'sequence', RBI_5UTR, 'format', 'diagram'); %illlustrate the sequence diagram 

string = "The score for our RBI_5UTR versus the WT%d' is = %d";

%RBI_5UTR versus WT1
RBI_5UTR_v_WT1 = score_comp(RBI_5UTR_dot_bracket,wild_type_1);
WT_value_1 = 1; 
sprintf(string,WT_value_1,RBI_5UTR_v_WT1)

%RBI_5UTR versus WT2
RBI_5UTR_v_WT2 = score_comp(RBI_5UTR_dot_bracket,wild_type_2);
WT_value_2 = 2; 
sprintf(string,WT_value_2,RBI_5UTR_v_WT2)

%RBI_5UTR versus WT3
RBI_5UTR_v_WT3 = score_comp(RBI_5UTR_dot_bracket,wild_type_3);
WT_value_3 = 3; 
sprintf(string,WT_value_3,RBI_5UTR_v_WT3)




%step 3 

RBI_5UTR_table = step_3(RBI_5UTR); %called our function with our RBI_5UTR string 
 
writetable(RBI_5UTR_table,'RB1_mdubois6.txt'); %wrote it to a .txt file
type RB1_mdubois6.txt;
 
[min_value, indx_value] = (min(RBI_5UTR_table{:,4})); %found the min value and index of the min value
 
step4_position = RBI_5UTR_table{indx_value,1}; %used the index to find the position 
step4_og_nucleotide = RBI_5UTR_table{indx_value,2}; %used the index to find the og nucleotide 
step4_mutated_nucleotide = RBI_5UTR_table{indx_value,3}; %used the index to find the mutated nucleotide

RBI_5UTR_for_graph = RBI_5UTR; %created a new sequence to hold this single mutation that I wanted to plot 
RBI_5UTR_for_graph(step4_position) = step4_mutated_nucleotide; %executed this single mutation 
 
mutated_structure = rnaplot(rnafold(RBI_5UTR_for_graph, 'MinLoopSize',7),'sequence', RBI_5UTR_for_graph, 'format', 'diagram'); %plotted the single mutation 

step_4 = 'The most different strucutre of RBI_5UTR had a comparison score of %d at position %d where %s was mutated into %s';
step_4_to_display = sprintf(step_4,min_value,step4_position,step4_og_nucleotide,step4_mutated_nucleotide)





%step 5 sequence 

step_5_dna = upper('acgcaaaagaaggcaagatctcttttttcttttgtgttgtcatatac'); %my DNA sequence converted to all caps
step_5_rna = dna2rna(step_5_dna) %converted to RNA

step_5_dot_bracket = rnafold(step_5_rna, 'MinLoopSize',7); %turned into dot bracket
step_5_circle = rnaplot(step_5_dot_bracket, 'sequence', step_5_rna, 'format', 'diagram'); %turned into its diagram

hold off

step_5_table = step_3(step_5_rna); % ran my function against my new sequence 

[step_5_min_value, step_5_indx_value] = (min(step_5_table{:,4})); %repeated all the same steps as 
 
step5_position = step_5_table{step_5_indx_value,1};
step5_og_nucleotide = step_5_table{step_5_indx_value,2};
step5_mutated_nucleotide = step_5_table{step_5_indx_value,3};

step5_for_graph = RBI_5UTR;
step5_for_graph(step4_position) = step4_mutated_nucleotide;

step_5 = 'The most different strucutre of rs4141463 had a comparison score of %d at position %d where %s was mutated into %s';
step_5_to_display = sprintf(step_5,step_5_min_value,step5_position,step5_og_nucleotide,step5_mutated_nucleotide)


writetable(step_5_table,'step5_mdubois6.txt');
type step5_mdubois6.txt;

first_codon_locations = [1,4,7,10,13,16,19,22,25,28,31,34,37,40,43]; %found all the positions of the first/second/third codons, yes the long way 
second_codon_locations = [2,5,8,11,14,17,20,23,26,29,32,35,38,41,44];
third_codon_locations = [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45];

first_codon_cscores = step_5_table{:,4}(ismember(step_5_table{:,1},first_codon_locations));
second_codon_cscores = step_5_table{:,4}(ismember(step_5_table{:,1},second_codon_locations));
third_codon_cscores = step_5_table{:,4}(ismember(step_5_table{:,1},third_codon_locations));

ave_first_codon_comparison_scores = mean(first_codon_cscores); %ismember allowed me to create a boolean selector to select for all the comparison scores at the 1st codon position  
ave_second_codon_comparison_scores = mean(second_codon_cscores); %repeated the above for the second and third codon position 
ave_third_codon_comparison_scores = mean(third_codon_cscores); 

max_first_codon_comparison_scores = max(first_codon_cscores); %ismember allowed me to create a boolean selector to select for all the comparison scores at the 1st codon position  
max_second_codon_comparison_scores = max(second_codon_cscores); %repeated the above for the second and third codon position 
max_third_codon_comparison_scores = max(third_codon_cscores); 

min_first_codon_comparison_scores = min(first_codon_cscores); %ismember allowed me to create a boolean selector to select for all the comparison scores at the 1st codon position  
min_second_codon_comparison_scores = min(second_codon_cscores); %repeated the above for the second and third codon position 
min_third_codon_comparison_scores = min(third_codon_cscores); 

std_1 = std(first_codon_cscores);
std_2 = std(second_codon_cscores);
std_3 = std(third_codon_cscores);

step_6_string = ('Average comparison score for mutations in the 1st codon = %d, 2nd codon = %d, 3rd codon = %d');
sprintf(step_6_string,ave_first_codon_comparison_scores,ave_second_codon_comparison_scores,ave_third_codon_comparison_scores)


function [T] = step_3(sequence)
    og_sequence = sequence; 
    sequence_dot_bracket = rnafold(sequence, 'MinLoopSize',7); %Step 1: Illustrate the secondary structure of your sequence 

    T = table('Size',[0 4],'VariableTypes',{'single','string','string','single'},'VariableNames',{'Mutated Position','Original Nucleotide at this position','Nucleotide was mutated to this nucleotide','Comparison score (higher is more similar)'}); %initiated a table with no rows....yet, 4 columns all labeled accordingly with their column titles as well 
    j = 1; %my counter j which indicates which row we are on 

    for i = 1:length(og_sequence) %looped through the entire length of the sequence
        func_sequence = sequence; % resetting my sequence so that once i am done with a certain nucleotide it goes back to being the og

        if func_sequence(i) == 'G' %run through every possible combination of nucleotides 
            func_sequence(i) = 'C'; %mutate the nucleotide  
            dot_bracket = rnafold(func_sequence); %fold the mutated sequence 
            score = score_comp(sequence_dot_bracket,dot_bracket); %calculate the comparison score between the mutated sequence and the og sequence 
            T(j,:) = {i,'G','C',score}; %add all this info to the particular row in the table 

            j = j+1; %add 1 to j so that we know to go to the next row and repeat the process 
            func_sequence(i) = 'U'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'G','U',score};

            j = j+1;
            func_sequence(i) = 'A'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'G','A',score};
            j = j+1;

        elseif func_sequence(i) == 'C'
            func_sequence(i) = 'G'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'C','G',score};

            j = j+1;
            func_sequence(i) = 'U'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'C','U',score};

            j = j+1;
            func_sequence(i) = 'A'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'C','A',score};
            j = j+1; 

        elseif func_sequence(i) == 'A'
            func_sequence(i) = 'G'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'A','G',score};

            j = j+1;
            func_sequence(i) = 'U'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'A','U',score};

            j = j+1;
            func_sequence(i) = 'C'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'A','C',score};
            j = j+1;

        elseif func_sequence(i) == 'U'
            func_sequence(i) = 'G'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'U','G',score};

            j = j+1;
            func_sequence(i) = 'C'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'U','C',score};

            j = j+1;
            func_sequence(i) = 'A'; 
            dot_bracket = rnafold(func_sequence);
            score = score_comp(sequence_dot_bracket,dot_bracket);
            T(j,:) = {i,'U','A',score};
            j = j+1;
        end 
    end
end


function [ score ] = score_comp( str1, str2 )
%SCORE_COMP compares two RNA structures in dot-bracket notation
%   Takes two dot-bracket notation RNA secondary structures as input and
%   compares their structure by aligning them and outputting a score as a
%   function of the RNA's length. The highest score of 1 indicates the same
%   structure and lower scores (lowest is 0) indicate dissimilar
%   structures.

str1 = strrep(str1,'.','A'); str1 = strrep(str1,'(','C'); str1 = strrep(str1,')','G');
str2 = strrep(str2,'.','A'); str2 = strrep(str2,'(','C'); str2 = strrep(str2,')','G');

score = nwalign(str1,str2,'Alphabet','NT','ScoringMatrix',eye(4),'GAPOPEN',.5,'EXTENDGAP',0.25);
score = score / max(length(str1),length(str2));

end


