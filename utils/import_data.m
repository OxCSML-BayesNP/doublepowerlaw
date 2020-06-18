function [wordcounts,words] = import_data(dataset)

rep_data = '../data/'; % folder of the data

switch dataset
    case 'ANC'
        namedataset = 'American National Corpus';
        matfilename = [rep_data 'ANC/ANC.mat'];
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            % To create mat file
            filename = './ANC/ANC-written-count.txt';
            A = importdata([rep_data filename]);
            words = A.textdata(2:end, 1); % list of words
            wordcounts = A.data; % vector of wordcounts
            save(matfilename, 'wordcounts', 'words')
        end        
    case 'googlengram'
        namedataset = 'Googlengram';
        matfilename = [rep_data 'googlengram/googlengram.mat']; 
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            % To create mat file
            filename = './googlengram/frequency-alpha-oed.xlsx';        
            [Anum, Atext] = xlsread([rep_data filename]);        
            words = Atext(2:end, 2);
            wordcounts = Anum(:, 3);
            matfilename = [rep_data 'googlengram/googlengram.mat']; 
            save(matfilename, 'wordcounts', 'words')       
        end
    case 'nipsold'
        namedataset = 'nipsold';
        matfilename =  [rep_data 'nipsold/nips.mat']; 
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            A=importdata('../data/nips/NIPS_1987-2015.csv');
            words = A.textdata(2:end,1);
            wordpaper = A.data;
            wordcounts = sum(wordpaper, 2);
            [~, ind] = sort(wordcounts, 'descend');
            words = words(ind);
            wordcounts = wordcounts(ind);
            wordpaper = wordpaper(ind, :);
            save(matfilename, 'wordcounts', 'words', 'wordpaper')
        end
    case 'nips'
        namedataset = 'nips';
        matfilename =  [rep_data 'nips/nips.mat']; 
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            load([rep_data 'nips/nipspapersmatrix.mtx'])
            wordcounts = spconvert(nipspapersmatrix(2:end, :));

            words = importdata([rep_data 'nips/nipspaperswords.txt']);
            papers = importdata([rep_data 'nips/nipspapersarticles.txt']);
            for i=1:length(words)
                words{i}=words{i}(2:end-1);
            end
            for i=1:length(papers)
                papers{i}=papers{i}(2:end-1);
            end
            save(matfilename, 'wordcounts', 'words', 'papers')
        end
        wordcounts = sum(wordcounts, 1);
    case 'nipsNSW'
        namedataset = 'nipsNSW';
        matfilename =  [rep_data 'nips/nipsNSW.mat']; 
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            load([rep_data 'nips/nipspapersmatrixNSW.mtx'])
            wordcounts = spconvert(nipspapersmatrixNSW(2:end, :));

            words = importdata([rep_data 'nips/nipspaperswordsNSW.txt']);
            papers = importdata([rep_data 'nips/nipspapersarticlesNSW.txt']);
            for i=1:length(words)
                words{i}=words{i}(2:end-1);
            end
            for i=1:length(papers)
                papers{i}=papers{i}(2:end-1);
            end
            save(matfilename, 'wordcounts', 'words', 'papers')
        end
        wordcounts = sum(wordcounts, 1); 
    case 'nips1000'
        load([rep_data 'nips/nips1000.mat']);
        wordcounts = cnts';
        words = [];
    case 'superusers'
        wordcounts = dlmread([rep_data 'Superusers/superusers_cnts.txt']);
        wordcounts = wordcounts';
        words=[];
    case 'travian'
        wordcounts = dlmread([rep_data 'Travian/travian_mess_1_13.txt']);
        wordcounts = wordcounts';
        words=[];
    case 'simulated_m1'
        load([rep_data 'Simulated_m1/simulated_data.mat']);
        wordcounts = cnts';
        words=[];
    case 'englishbooks'
        namedataset = 'englishbooks';
        n_books = 20;
        matfilename =  [rep_data 'englishbooks/englishbooks.mat'];
        load(matfilename);
        wordcounts = sum(data.frequencies(1:n_books,:), 1)';
        idx_words = find(wordcounts>0);
        wordcounts = wordcounts(idx_words)';
        words = data.words(idx_words);
        % sort
        [~, ind] = sort(wordcounts, 'descend');
            words = words(ind);
            wordcounts = wordcounts(ind)';
    case 'trumptweets'
        namedataset = 'trumptweets';
        matfilename =  [rep_data 'trumptweets/trumptweets.mat']; 
        if exist(matfilename, 'file') == 2
            load(matfilename);
        else
            load([rep_data 'trumptweets/trumptweetwordmatrix.mtx'])
            wordcounts = spconvert(trumptweetwordmatrix(2:end, :));
            words = importdata([rep_data 'trumptweets/trumptweetwords.txt']);
            for i=1:length(words)
                words{i}=words{i}(2:end-1);
            end
            save(matfilename, 'wordcounts', 'words')
        end
        wordcounts = sum(wordcounts, 1);
    otherwise
        error('Unknown dataset');
end

