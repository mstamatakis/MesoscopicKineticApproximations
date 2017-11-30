function out1 = Hamiltonian(varargin)

ApproxIdent = varargin{1};
H0 = varargin{2};
H1 = varargin{3};
h = varargin{4};
nn = varargin{5};

switch ApproxIdent
    case 'BP'
        g1 = varargin{6};
        epsl = varargin{7};
        out1 = HamiltonianBP(H0,H1,h,nn,g1,epsl);
    case 'BPE'
        g1 = varargin{6};
        epsl = varargin{7};
        out1 = HamiltonianBPE(H0,H1,h,nn,g1,epsl);
    case 'BPEC'
        g1 = varargin{6};
        c1 = varargin{7};
        epsl = varargin{8};
        out1 = HamiltonianBPEC(H0,H1,h,nn,g1,c1,epsl);
    case 'K2NNC1'
        g1 = varargin{6};
        g2 = varargin{7};
        c1 = varargin{8};
        c2 = varargin{9};
        epsl = varargin{10};
        out1 = HamiltonianK2NNC1(H0,H1,h,nn,g1,g2,c1,c2,epsl);
    case 'K2NNC2'
        g1 = varargin{6};
        g2 = varargin{7};
        c1 = varargin{8};
        c2 = varargin{9};
        p1 = varargin{10};
        p2 = varargin{11};
        epsl = varargin{12};
        out1 = HamiltonianK2NNC2(H0,H1,h,nn,g1,g2,c1,c2,p1,p2,epsl);
    case 'K2NNC2T1'
        g1 = varargin{6};
        g2 = varargin{7};
        c1 = varargin{8};
        c2 = varargin{9};
        p1 = varargin{10};
        p2 = varargin{11};
        p3 = varargin{12};
        epsl = varargin{13};
        out1 = HamiltonianK2NNC2T1(H0,H1,h,nn,g1,g2,c1,c2,p1,p2,p3,epsl);
    case 'K2NNC3'
        g1 = varargin{6};
        g2 = varargin{7};
        c1 = varargin{8};
        c2 = varargin{9};
        p1 = varargin{10};
        p2 = varargin{11};
        p3 = varargin{12};
        epsl = varargin{13};
        out1 = HamiltonianK2NNC3(H0,H1,h,nn,g1,g2,c1,c2,p1,p2,p3,epsl);
    case 'K3NNC1'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        epsl = varargin{13};
        out1 = HamiltonianK3NNC1(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,epsl);      
    case 'K3NN'
        g1 = varargin{6};
        g2 = varargin{7};
        c1 = varargin{8};
        c2 = varargin{9};
        c3 = varargin{10};
        epsl = varargin{11};
        out1 = HamiltonianK3NN(H0,H1,h,nn,g1,g2,c1,c2,c3,epsl);
    case 'K3NN1corr'
        g1 = varargin{6};
        epsl = varargin{7};
        out1 = HamiltonianK3NN1corr(H0,H1,h,nn,g1,epsl);
    case 'K3NN2corr'
        g1 = varargin{6};
        g2 = varargin{7};
        epsl = varargin{8};
        out1 = HamiltonianK3NN2corr(H0,H1,h,nn,g1,g2,epsl);
    case 'K3NNC2'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        p1 = varargin{13};
        p2 = varargin{14};
        p3 = varargin{15};
        epsl = varargin{16};
        out1 = HamiltonianK3NNC2(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,epsl);
    case 'K3NNC3'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        p1 = varargin{13};
        p2 = varargin{14};
        p3 = varargin{15};
        p4 = varargin{16};
        p5 = varargin{17};
        p6 = varargin{18};
        epsl = varargin{19};
        out1 = HamiltonianK3NNC3(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,epsl);
    case 'K3NNC3b'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        p1 = varargin{13};
        p2 = varargin{14};
        p3 = varargin{15};
        p4 = varargin{16};
        p5 = varargin{17};
        p6 = varargin{18};
        epsl = varargin{19};
        out1 = HamiltonianK3NNC3b(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,epsl);
    case 'K3NNC4'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        p1 = varargin{13};
        p2 = varargin{14};
        p3 = varargin{15};
        p4 = varargin{16};
        p5 = varargin{17};
        p6 = varargin{18};
        p7 = varargin{19};
        p8 = varargin{20};
        epsl = varargin{21};
        out1 = HamiltonianK3NNC4(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,p7,p8,epsl);
    case 'K3NNC5'
        g1 = varargin{6};
        g2 = varargin{7};
        g3 = varargin{8};
        c1 = varargin{9};
        c2 = varargin{10};
        c3 = varargin{11};
        c4 = varargin{12};
        p1 = varargin{13};
        p2 = varargin{14};
        p3 = varargin{15};
        p4 = varargin{16};
        p5 = varargin{17};
        p6 = varargin{18};
        p7 = varargin{19};
        p8 = varargin{20};
        p9 = varargin{21};
        epsl = varargin{22};
        out1 = HamiltonianK3NNC5(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,p7,p8,p9,epsl);
    otherwise
        error('Unknown approximate Hamiltonian!')
end

end