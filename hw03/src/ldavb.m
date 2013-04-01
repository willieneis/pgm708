function out = ldavb(data,b,alpha)
% lda variational bayes inference, given data matrix, pre-inferred
%   topics b, and alpha dir param.

    numTopics = size(b,2);
    numDocs = size(data,1);
    numVocab = size(data,2);
    gam_init = ones(1,numTopics);
    iters = zeros(1,numDocs);
    epsval = 0.001;

    % do inference for each doc (row in data)
    for i=1:numDocs
        gam = gam_init;
        phi = zeros(sum(data(i,:)),numTopics);
        docVec = convertDoc(data(i,:));
        eps_gam = Inf; eps_phi = Inf;
        while eps_gam>epsval || eps_phi>epsval
            phi_old = phi;
            gam_old = gam;
            phi = update_phi(size(phi,1));
            gam = update_gam();
            eps_phi = max(max(abs(phi-phi_old)));
            eps_gam = max(abs(gam-gam_old));
            iters(i) = iters(1)+1;
        end
        gamMat(i,:) = gam;
        phiCell{i} = phi;
        fprintf('Finished doc %d.\n',i);
    end
    out = {gamMat,phiCell,iters};

    % -----------

    function convDoc = convertDoc(doc)
        convDoc = [];
        for n=1:length(doc)
            if doc(n)>0
                convDoc(end+1:end+doc(n)) = n;
            end
        end
    end

    function gam = update_gam()
        gam = alpha+sum(phi);
    end

    function phi = update_phi(numWords)
        for w=1:numWords
            obsWord = docVec(w);
            topicRow = b(obsWord,:);
            phi(w,:) = topicRow.*exp(psi(gam)-psi(sum(gam)));
            phi(w,:) = phi(w,:)/sum(phi(w,:));
        end
    end

end
