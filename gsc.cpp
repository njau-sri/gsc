#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "vcfio.h"

using std::size_t;

using bint = int;

extern "C"
{
    void dsyevr_(char *jobz, char *range, char *uplo, bint *n, double *a, bint *lda, double *vl, double *vu,
        bint *il, bint *iu, double *abstol, bint *m, double *w, double *z, bint *ldz, bint *isuppz,
        double *work, bint *lwork, bint *iwork, bint *liwork, bint *info);

} // extern "C"

namespace
{
    struct Par
    {
        std::string vcf;
        std::string out;
        int top = 10;
    };

    Par par;

    template<typename T>
    int count_shared(T a, T b, T c, T d)
    {
        if (a == b)
            return static_cast<int>(a == c) + (a == d);

        if (c == d)
            return static_cast<int>(c == a) + (c == b);

        return static_cast<int>(a == c) + (a == d) + (b == c) + (b == d);
    }

    void calc_gsc_1(int n, const std::vector<allele_t> &x, const std::vector<allele_t> &y, int &a, int &b)
    {
        for (int i = 0; i < n; ++i) {
            if (x[i] != 0 && y[i] != 0) {
                ++a;
                if (x[i] == y[i])
                    ++b;
            }
        }
    }

    void calc_gsc_2(int n, const std::vector<allele_t> &x, const std::vector<allele_t> &y, int &a, int &b)
    {
        for (int i = 0; i < n; ++i) {
            int j = i * 2, k = j + 1;
            if (x[j] != 0 && x[k] != 0 && y[j] != 0 && y[k] != 0) {
                a += 2;
                b += count_shared(x[j], x[k], y[j], y[k]);
            }
        }
    }

    void calc_gsc_matrix(const Genotype &gt, std::vector< std::vector<double> > &x)
    {
        static const int nb = 10000;

        int m = 0;
        int n = gt.ind.size();
        std::vector< std::vector<int> > z(n, std::vector<int>(n, 0));

        x.assign(n, std::vector<double>(n, 0.0));

        if (gt.ploidy == 1) {
            std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(nb));

            for (auto &v : gt.dat) {
                if (m < nb) {
                    for (int i = 0; i < n; ++i)
                        dat[i][m] = v[i];
                    ++m;
                    continue;
                }
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j)
                        calc_gsc_1(m, dat[i], dat[j], z[i][j], z[j][i]);
                for (int i = 0; i < n; ++i)
                    dat[i][0] = v[i];
                m = 1;
            }

            if (m > 0) {
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j)
                        calc_gsc_1(m, dat[i], dat[j], z[i][j], z[j][i]);
            }
        }
        else {
            std::vector< std::vector<allele_t> > dat(n, std::vector<allele_t>(nb * 2));

            for (auto &v : gt.dat) {
                if (m < nb) {
                    for (int i = 0; i < n; ++i) {
                        dat[i][m * 2] = v[i * 2];
                        dat[i][m * 2 + 1] = v[i * 2 + 1];
                    }
                    ++m;
                    continue;
                }
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j)
                        calc_gsc_2(m, dat[i], dat[j], z[i][j], z[j][i]);
                for (int i = 0; i < n; ++i) {
                    dat[i][0] = v[i * 2];
                    dat[i][1] = v[i * 2 + 1];
                }
                m = 1;
            }

            if (m > 0) {
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j)
                        calc_gsc_2(m, dat[i], dat[j], z[i][j], z[j][i]);
            }
        }

        for (int i = 0; i < n; ++i) {
            x[i][i] = 1.0;
            for (int j = i + 1; j < n; ++j) {
                if (z[i][j] != 0) {
                    auto a = static_cast<double>(z[j][i]) / z[i][j];
                    x[i][j] = x[j][i] = a;
                }
            }
        }
    }

    int syevr(char jobz, char range, char uplo, bint n, double *a, bint lda, double vl, double vu, bint il, bint iu,
        double abstol, bint *m, double *w, double *z, bint ldz, bint *isuppz)
    {
        bint info = 0;

        double wkopt;
        bint lwork = -1, iwkopt, liwork = -1;
        dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
            &wkopt, &lwork, &iwkopt, &liwork, &info);

        if (info == 0) {
            lwork = static_cast<bint>(wkopt);
            liwork = iwkopt;
            bint *iwork = new bint[liwork];
            double *work = new double[lwork];
            dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz,
                work, &lwork, iwork, &liwork, &info);
            delete[] work;
            delete[] iwork;
        }

        return info;
    }

    void eigen(const std::vector< std::vector<double> > &mat, std::vector<double> &eval, std::vector< std::vector<double> > &evec)
    {
        std::vector<double> a;
        for (auto &v : mat)
            a.insert(a.end(), v.begin(), v.end());

        int n = mat.size();
        std::vector<double> w(n), z(n*n);
        std::vector<bint> sup(2 * n);
        bint m = 0;

        syevr('V', 'A', 'U', n, a.data(), n, 0.0, 0.0, 0, 0, 0.0, &m, w.data(), z.data(), n, sup.data());

        std::reverse(w.begin(), w.end());

        eval.swap(w);
        evec.clear();

        for (int j = 0; j < n; ++j) {
            auto k = n - 1 - j;
            std::vector<double> v(n);
            for (int i = 0; i < n; ++i)
                v[i] = z[k*n + i];
            evec.push_back(v);
        }
    }

    int write_gsc_matrix(const std::vector<std::string> &ind, const std::vector< std::vector<double> > &mat, const std::string &filename)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
            return 1;
        }

        int n = ind.size();

        for (int i = 0; i < n; ++i) {
            ofs << ind[i];
            for (int j = 0; j < n; ++j)
                ofs << "\t" << mat[i][j];
            ofs << "\n";
        }

        return 0;
    }

    int write_eigenvalues(std::vector<double> &eval, const std::string &filename)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
            return 1;
        }

        double sum = std::accumulate(eval.begin(), eval.end(), 0.0);

        for (auto e : eval)
            ofs << e << "\t" << 100 * e / sum << "\n";

        return 0;
    }

    int write_eigenvectors(const std::vector<std::string> &ind, const std::vector< std::vector<double> > &evec, const std::string &filename)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
            return 1;
        }

        int n = ind.size();
        int m = evec.size();

        ofs << "Indiv";
        for (int j = 0; j < m; ++j)
            ofs << "\tEV" << j + 1;
        ofs << "\n";

        for (int i = 0; i < n; ++i) {
            ofs << ind[i];
            for (int j = 0; j < m; ++j)
                ofs << "\t" << evec[j][i];
            ofs << "\n";
        }

        return 0;
    }

} // namespace

int gsc(int argc, char *argv[])
{
    std::cerr << "gsc (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd("gsc [options]");

    cmd.add("--vcf", "VCF file", "");
    cmd.add("--out", "output file", "gsc.out");
    cmd.add("--top", "number of eigenvectors", "10");

    if (argc < 2) {
        cmd.help();
        return 1;
    }

    cmd.parse(argc, argv);

    par.vcf = cmd.get("--vcf");
    par.out = cmd.get("--out");
    par.top = std::stoi(cmd.get("--top"));

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::vector< std::vector<double> > mat;

    calc_gsc_matrix(gt, mat);

    std::vector<double> eval;
    std::vector< std::vector<double> > evec;

    eigen(mat, eval, evec);

    eval.resize(par.top);
    evec.resize(par.top);

    write_gsc_matrix(gt.ind, mat, par.out + ".mat");

    write_eigenvalues(eval, par.out + ".eval");

    write_eigenvectors(gt.ind, evec, par.out + ".evec");

    return 0;
}
