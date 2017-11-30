function out1 = AllStates(nsites)

persistent AllStates1 AllStates7 AllStates13 AllStates19

switch nsites
    case 1
        computeanew = size(AllStates1,1) == 0;
    case 7
        computeanew = size(AllStates7,1) == 0;
    case 13
        computeanew = size(AllStates13,1) == 0;
    case 19
        computeanew = size(AllStates19,1) == 0;
    otherwise
        error(['AllStates function called for an invalid ' ...
            'number of sites = ' num2str(nsites) ...
            '. nsites must be one of ' ...
            '1, 7, 13, 19.']);
end
        

if computeanew
    
    disp(['Computing AllStates matrix for nsites = ' num2str(nsites) '...']);

    AllStatesTemp = de2bi(0:2^nsites-1,nsites); 
    AllStatesTemp = AllStatesTemp(:,end:-1:1);

    disp('AllStates matrix created!');

    eval(['AllStates' num2str(nsites) ' = AllStatesTemp;']);
end

eval(['out1 = AllStates' num2str(nsites) ';']);

return