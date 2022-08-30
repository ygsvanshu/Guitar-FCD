#include <iostream>
#include <ctime>
#include "Eigen/Sparse"

int main (int argc, char *argv[]) {

    clock_t t0,t1;
    double ts;
    double dx,dy;

    int nx,ny,nn;
    int im,ic,ip;
    int jm,jc,jp;
    int row,col;
    int n,niter;

    nx = 341;
    ny = 341;
    nn = nx*ny;

    dx = 1.0; //1.18533E-3;
    dy = 1.0; //1.18533E-3;

    niter = 10;

    std::cout << std::endl;

    t0 = clock();

    Eigen::SparseMatrix<double> A0(2*nn,nn);
    Eigen::SparseMatrix<double> AT(nn,2*nn);
    Eigen::SparseMatrix<double> AA(nn,nn);

    std::cout << Eigen::Success << std::endl;
    std::cout << Eigen::NumericalIssue << std::endl;
    std::cout << Eigen::NoConvergence << std::endl;
    std::cout << Eigen::InvalidInput << std::endl;

    A0.reserve(Eigen::VectorXi::Constant(nn,4));

    row = 0;
    for (ic=0;ic<nx;ic++){
        im = std::max(0,ic-1);
        ip = std::min(nx-1,ic+1);
        for (jc=0;jc<ny;jc++){
            
            col = (jc*nx)+im;
            A0.insert(row,col) = 1.0/(double(im-ip)*dx);

            col = (jc*nx)+ip;
            A0.insert(row,col) = 1.0/(double(ip-im)*dx);

            row++;
        }
    }

    for (jc=0;jc<ny;jc++){
        jm = std::max(0,jc-1);
        jp = std::min(ny-1,jc+1);
        for (ic=0;ic<nx;ic++){
            
            col = (jm*nx)+ic;
            A0.insert(row,col) = 1.0/(double(jm-jp)*dy);

            col = (jp*nx)+ic;
            A0.insert(row,col) = 1.0/(double(jp-jm)*dy);
            
            row++;
        }
    }

    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"A0 built in "<<ts<<" seconds"<<std::endl;

    t0 = clock();
    A0.makeCompressed();
    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"A0 compressed in "<<ts<<" seconds"<<std::endl;

    t0 = clock();
    AT = A0.transpose();
    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"A0 transposed in "<<ts<<" seconds"<<std::endl;

    t0 = clock();
    AA = AT*A0;
    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"AA = AT*A0 computed in "<<ts<<" seconds"<<std::endl;

    t0 = clock();
    // Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    // Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
    Eigen::BiCGSTAB < Eigen::SparseMatrix<double> > solver;
    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"Solver built in "<<ts<<" seconds"<<std::endl;

    t0 = clock();
    // solver.compute(A0);
    solver.compute(AA);
    t1 = clock();
    ts = double(t1-t0)/CLOCKS_PER_SEC;
    std::cout<<"Solver factorized in "<<ts<<" seconds"<<std::endl;
    
    Eigen::VectorXd B0 = Eigen::VectorXd::Random(2*nn);
    Eigen::VectorXd BB(nn);

    Eigen::VectorXd X;

    // std::cout<<"Sparse matrix AA has "<<AA.nonZeros()<<" values "<<AA.isCompressed()<<std::endl;

    for (n=0;n<niter;n++){

        // for (row=0;row<(2*nn);row++){
        //     B(row) = double(rand()/RAND_MAX);
        // }

        // B.random();

        t0 = clock();
        BB = AT*B0;
        X = solver.solve(BB);
        t1 = clock();
        ts = double(t1-t0)/CLOCKS_PER_SEC;
        std::cout<<"Solved instance " << n << " Iterations = " << solver.iterations() << " Error = " << solver.error() << " Info = " << solver.info() << " ClockTime = "<<ts<<" seconds"<<std::endl;
        // std::cout<<"Solved instance " << n << " ClockTime = "<<ts<<" seconds"<<std::endl;
    }

    std::cout << std::endl;

}