import numpy as np
import scipy
import scipy.linalg
from scipy.special import eval_chebyt
from sklearn import linear_model
import itertools
import warnings


class ChIMES:
    def __init__(self):
        return
        
    def tersoff_smooth(self, crds, morse_fo, crds_out):
        y = np.zeros(crds.shape)
        dt = crds_out * (1 - morse_fo)

        mask_out = crds > crds_out
        mask_in = crds < dt
        mask_between = ~(mask_out | mask_in)

        y[mask_in] = 1.0

        frac = (crds[mask_between] - dt) / (crds_out - dt)
        y[mask_between] = 0.5 + 0.5 * np.sin(np.pi * frac + 0.5 * np.pi)
        return y

    def morse_trans(self, crds, crds_in, crds_out, morse_lambda):
        x = np.exp(-crds / morse_lambda)
        x_in = np.exp(-crds_in / morse_lambda)
        x_out = np.exp(-crds_out / morse_lambda)
        return x, x_in, x_out

    def inverse_trans(self, crds, crds_in, crds_out):
        x = 1/crds
        x_in = 1/crds_in
        x_out = 1/crds_out
        return x, x_in, x_out

    def rescale_into_s(self, x, x_in, x_out, out_coeff=1.0):
        x_avg = (x_in + x_out*out_coeff) / 2
        x_diff = np.abs(x_in - x_out*out_coeff) / 2
        s = (x - x_avg) / x_diff
        return s

    def make_Amatrix(self, s, O2b, smooth_f, N_particles):
        n_datapoints = s.shape[0]
        assert n_datapoints == smooth_f.shape[0], "number of data points does not match"

        A = []
        for o in range(1, O2b + 1):
            column = eval_chebyt(o, s) * smooth_f
            A.append(column)
        A.append(np.ones_like(s) * N_particles)

        return np.column_stack(A)
        
    def make_Amatrix_design(self, s, O2b, s_lamb, Olamb, smooth_fr, N_particles):
        n_datapoints = s.shape[0]
        assert n_datapoints == s_lamb.shape[0], "number of data points does not match"
        assert n_datapoints == smooth_fr.shape[0], "number of data points does not match"

        A = []
        for o_lamb in range(0, Olamb):
            for o_r in range(1, O2b+1):
                column = eval_chebyt(o_r, s) * eval_chebyt(o_lamb, s_lamb) * smooth_fr
                A.append(column)
        A.append(np.ones_like(s) * N_particles)

        return np.column_stack(A)

    def make_3bAmatrix(
            self,
            s_ij,
            s_ik,
            s_jk,
            smooth_f_ij,
            smooth_f_ik,
            smooth_f_jk,
            O2b,
            O3b,
            N_particles
    ):
        """
        Note: eval_chebyt can still return extrapolation value when s is outside
        the interval of [-1, 1]. However, the smoothing function guarantee that 
        the extrapolation value is zero out so it doesn't influence the A matrix
        generation.
        """
        n_datapoints = s_ij.shape[0]
        assert s_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert s_jk.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_f_ij.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_f_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_f_jk.shape[0] == n_datapoints, "number of data points does not match"

        print_out_list = []
        A = []

        print_out_list.append("Two-body order:")
        for o in range(1, O2b + 1):
            print_out_list.append(str(o))
            column = eval_chebyt(o, s_ij) * smooth_f_ij + eval_chebyt(o, s_jk) * smooth_f_jk + eval_chebyt(o,
                                                                                                           s_ik) * smooth_f_ik
            A.append(column)

        print_out_list.append("Three-body order & equivalent terms:")
        count_string_list = [str(i) for i in range(0, O3b)]
        desire_order = ''.join(count_string_list)

        combinations = list(itertools.combinations_with_replacement(desire_order, 3))[O3b:]
        
        for irr_combination in combinations:
            possible_terms = list(set(
                list(itertools.permutations(irr_combination[0] + irr_combination[1] + irr_combination[2], 3))))
            column = 0
            print_out_list.append(str(possible_terms[0]) + " " + str(len(possible_terms)))
            for term in possible_terms:
                n1 = float(term[0])
                n2 = float(term[1])
                n3 = float(term[2])
                column += smooth_f_ij * smooth_f_ik * smooth_f_jk * eval_chebyt(n1, s_ij) * eval_chebyt(n2,
                                                                                                        s_ik) * eval_chebyt(
                    n3, s_jk)
            A.append(column)

        A.append(np.ones_like(s_ij) * N_particles)

        for item in print_out_list:
            print(item)
        return np.column_stack(A)

    def make_3bAmatrix_design(
            self,
            s_ij,
            s_ik,
            s_jk,
            lamb,
            smooth_fr_ij,
            smooth_fr_ik,
            smooth_fr_jk,
            smooth_fl,
            O2b,
            O3b,
            Olamb_2b,
            Olamb_3b,
            N_particles
    ):
        """
        Note: eval_chebyt can still return extrapolation value when s is outside
        the interval of [-1, 1]. However, the smoothing function guarantee that 
        the extrapolation value is zero out so it doesn't influence the A matrix
        generation.
        """
        n_datapoints = s_ij.shape[0]
        assert s_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert s_jk.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_ij.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_jk.shape[0] == n_datapoints, "number of data points does not match"
        assert lamb.shape[0] == n_datapoints, "number of data points does not match"

        A = []
        for o_r in range(1, O2b + 1):
            for o_lamb in range(0, Olamb_2b):
                column = (eval_chebyt(o_r, s_ij) * eval_chebyt(o_lamb, lamb) * smooth_fr_ij +
                          eval_chebyt(o_r, s_jk) * eval_chebyt(o_lamb, lamb) * smooth_fr_jk +
                          eval_chebyt(o_r, s_ik) * eval_chebyt(o_lamb, lamb) * smooth_fr_ik)
                A.append(column)

        count_string_list = [str(i) for i in range(0, O3b)]
        desire_order = ''.join(count_string_list)

        combinations = list(itertools.combinations_with_replacement(desire_order, 3))[O3b:]
        for o_lamb in range(0, Olamb_3b):
            for irr_combination in combinations:
                possible_terms = list(set(
                    list(itertools.permutations(irr_combination[0] + irr_combination[1] + irr_combination[2], 3))))
                column = 0
                for term in possible_terms:
                    n1 = float(term[0])
                    n2 = float(term[1])
                    n3 = float(term[2])
                    column += smooth_fr_ij * smooth_fr_ik * smooth_fr_jk * eval_chebyt(n1, s_ij) * eval_chebyt(n2,
                                                                                                            s_ik) * eval_chebyt(
                        n3, s_jk) * eval_chebyt(o_lamb, lamb)
                A.append(column)

        A.append(np.ones_like(s_ij) * N_particles)

        return np.column_stack(A)

    def make_3bOnlyAmatrix_design(
            self,
            s_ij,
            s_ik,
            s_jk,
            lamb,
            smooth_fr_ij,
            smooth_fr_ik,
            smooth_fr_jk,
            O3b,
            Olamb_3b,
    ):
        """
        Note: eval_chebyt can still return extrapolation value when s is outside
        the interval of [-1, 1]. However, the smoothing function guarantee that 
        the extrapolation value is zero out so it doesn't influence the A matrix
        generation.
        """
        n_datapoints = s_ij.shape[0]
        assert s_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert s_jk.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_ij.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_ik.shape[0] == n_datapoints, "number of data points does not match"
        assert smooth_fr_jk.shape[0] == n_datapoints, "number of data points does not match"
        assert lamb.shape[0] == n_datapoints, "number of data points does not match"

        A = []
        count_string_list = [str(i) for i in range(0, O3b)]
        desire_order = ''.join(count_string_list)

        combinations = list(itertools.combinations_with_replacement(desire_order, 3))[O3b:]
        for o_lamb in range(0, Olamb_3b):
            for irr_combination in combinations:
                possible_terms = list(set(
                    list(itertools.permutations(irr_combination[0] + irr_combination[1] + irr_combination[2], 3))))
                column = 0
                for term in possible_terms:
                    n1 = float(term[0])
                    n2 = float(term[1])
                    n3 = float(term[2])
                    column += (smooth_fr_ij * smooth_fr_ik * smooth_fr_jk * eval_chebyt(n1, s_ij) *
                               eval_chebyt(n2, s_ik) * eval_chebyt(n3, s_jk) * eval_chebyt(o_lamb, lamb))
                A.append(column)

        return np.column_stack(A)
        
    def solve_LSQ_SVD(
        self, 
        A, 
        b, 
        svd_regularization_ratio=1e-12, 
        normal_eq=False,
        if_return_svd_results=False
        ):
        """
        Solve ordinary least square problem through TSVD

        Args:
            A: Configuration traning data with dimensions (n_configs, n_polynomials).
            b: Energy labeling data with dimensions (n_configs,).
            svd_regularization_ratio: TSVD regularization strength. Drop the pricipal 
                componenets and sigular values if the corresponfing singular values 
                samller than the maximum singular value * svd_regularization_ratio.
            if_return_svd_results: If retrun left and right singular vectors. 
                Defaults to False.
        
        Returns:
            c: Polynomial coefficients with dimension (n_polynomials,).

        """
        assert A.shape[0] == b.shape[0], "A and b must have the same row dimension (data points number)"
            
        if normal_eq:
            ATA = A.T @ A
            ATb = A.T @ b 
        
            U, sigma, VT = np.linalg.svd(ATA, full_matrices=False)
            drop_idx = sigma / sigma[0] < svd_regularization_ratio * svd_regularization_ratio
            inv_sigma = 1 / sigma
            inv_sigma[drop_idx] = 0.0
            c = VT.T @ np.diag(inv_sigma) @ U.T @ ATb
        else:
            U, sigma, VT = np.linalg.svd(A, full_matrices=False)
            drop_idx = sigma / sigma[0] < svd_regularization_ratio
            inv_sigma = 1 / sigma
            inv_sigma[drop_idx] = 0.0
            c = VT.T @ np.diag(inv_sigma) @ U.T @ b
        
        if if_return_svd_results:
            return c, U, sigma, VT
        else:
            return c

    def solve_L2_LSQ(self, A, b, gamma, mode="qr"):
        """
        Solve L2 regularized least square problem through QR factorization or Cholesky
        and then apply forward and back substitution.
        Default to use QR factorization.

        Args:
            A: Configuration traning data with dimensions (n_configs, n_polynomials).
            b: Energy labeling data with dimensions (n_configs,).
            alpha: L2 regularization strength.
            mode: Valid modes are "qr", "cholesky", "svd", "sklearn". Defaults to qr.
        
        Returns:
            c: Polynomial coefficients with dimension (n_polynomials,).
        """
        assert A.shape[0] == b.shape[0], "A and b must have the same row dimension (data points number)"
        
        # form normal equation to add regularization
        if mode != "svd" or "sklearn":
            ATA = A.T @ A
            ATb = A.T @ b
            ATA_reg = ATA
            ATA_reg.flat[:: A.shape[1] + 1] += gamma  # add gamma along the diagoanl elements
            
        if mode == "qr":
            # householder QR
            Q, R = np.linalg.qr(ATA_reg)

            # multiply Q.T on both side and perform back substitution
            c = scipy.linalg.solve_triangular(R, Q.T @ ATb)
        elif mode == "cholesky":
            # cholesky factorization
            L, low = scipy.linalg.cho_factor(ATA_reg)
            # forward and backward substitution
            c = scipy.linalg.cho_solve((L, low), ATb)
        elif mode == "svd":
            U, s, VT = np.linalg.svd(A, full_matrices=False)
            keep_mask = s > s.max() * 1e-15  # limit the maximum 2-norm condition number to be 1e-15, according to np.linalg.pinv
            s_truncated = s[keep_mask]
            rank = s_truncated.shape[0]

            d = s_truncated / (s_truncated*s_truncated + gamma)
            UT_b = U[:, :rank].T @ b
            d_UT_b = d * UT_b
            c = VT[:rank, :].T @ d_UT_b
        elif mode == "sklearn":
            reg = linear_model.Ridge(alpha=gamma)
            reg.fit(A, b)
            c = reg.coef_
        return c

    def solve_L1_LSQ_coordinate_descent(self, A, b, gamma):
        """
        Copy from the rk-lindsey/chimes_lsq Github repo
        """
        reg = linear_model.Lasso(alpha=gamma, fit_intercept=False, max_iter=100000)
        reg.fit(A, b)
        return reg.coef_

    def solve_L1_LSQ_LARS(self, A, b, gamma):
        """
        Copy from the rk-lindsey/chimes_lsq Github repo
        """
        reg = linear_model.LassoLars(
            alpha=gamma,
            fit_intercept=False,
            fit_path=False,
            verbose=True,
            max_iter=100000
        )
        reg.fit(A, b)
        return reg.coef_.ravel()

def parse_xyzf(mb_xyzf_fn, N_particles):                                                                                                
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pmf_data = np.genfromtxt(mb_xyzf_fn, skip_header=1, invalid_raise=False)[:, -1]
        particle_data = np.genfromtxt(mb_xyzf_fn, skip_header=2, invalid_raise=False)[:, 1:4]
       
    num_frame = int(particle_data.shape[0] / N_particles)
    particle_data = particle_data.reshape(num_frame, N_particles, 3)
    return pmf_data, particle_data
