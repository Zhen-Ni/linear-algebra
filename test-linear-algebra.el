;;; test-linear-algebra.el --- test linear algebra
;;; Commentary:
;;; This package tests linear-algebra.



;;; Code:

(load (concat default-directory "linear-algebra.el"))

(defun la--test-function ()
  "Private function."
  (let* ((v1 [1 2])
	 (v2 [3 4])
	 (v3 [5 6])
	 (m1 (vector v1 v2))
	 (m2 (vector v1 v2 v3))
	 (m2t [[1 3 5] [2 4 6]])
	 (m3 [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]])
	 (m4 [[5 6] [7 8]])
	 (m5 [[1 2 3] [4 5 6]])
	 (m6 [[3.0 2.0] [4.0 1.0]])
	 (m7 [[1 2 3] [4 5 6] [7 8 9]])
	 )
	
    (cl-assert (equal (la-at m1 0 1) 2) "error in la-at")
    (cl-assert (equal (la-at v3 0) 5) "error in la-at")
    (cl-assert (equal (la-shape 9) nil) t "error in la-shape")
    (cl-assert (equal (la-shape v1) '(2)) t "error in la-shape")
    (cl-assert (equal (la-shape m2) '(3 2)) t "error in la-shape")
    (cl-assert (equal (la-transpose m2) m2t) t "error in la-transpose")
    (cl-assert (equal (la-identity 4) m3) t "error in la-identity")
    (cl-assert (equal (la-diag m7) [1 5 9]) "error in la-diag")

    (cl-assert (equal (la-neg-v v1) [-1 -2]) "error in la-neg-v")
    (cl-assert (equal (la-neg-m m1) [[-1 -2] [-3 -4]]) "error in la-neg-m")
    (cl-assert (equal (la-cwise-vv * v1 [4 5]) [4 10]) "error in la-cwise-vv")
    (cl-assert (equal (la-cwise-mm + m1 m4) [[6 8] [10 12]]) "error in la-cwise-mm")
    (cl-assert (equal (la-mul-vv v1 v2) (la-cwise-vv * v1 v2)) "error in la-mul-vv")
    (cl-assert (equal (la-add-mm m1 m4) (la-cwise-mm + m1 m4)) "error in la-add-mm")
    (cl-assert (equal (la-dot-vv [1 2 3] [3 2 1]) 10) "error in la-dot-vv")
    (cl-assert (equal (la-dot-mv m1 [2 1]) [4 10]) "error in la-dot-vm")
    (cl-assert (equal (la-dot-vm [2 1] (la-transpose m1)) [4 10]) "error in la-dot-mv")
    (cl-assert (equal (la-dot-mm m5 [[1] [2] [3]]) [[14] [32]]) "error in la-dot-mm")
    (cl-assert (equal (la-kron [[1 2] [3 4]] [[1 2] [3 4]]) [[1 2 2 4] [3 4 6 8] [3 6 4 8] [9 12 12 16]]) "error in la-kron")

    (cl-assert (equal (la-norm1-m [[1 2] [3 -4]]) 10) "error in la-norm1-v")
    (cl-assert (equal (la-norm2-v [3 4]) 5.0) "error in la-norm2-v")
    (cl-assert (< (la-norm1-v (la-sub-vv (la-solve-mv m1 [4 10]) [2 1])) 1e-5) "error in la-solve-mv")
    (cl-assert (< (la-norm1-m (la-sub-mm (la-dot-mm m1 (la-inv m1)) (la-identity 2))) 1e-5) "error in la-solve-mv")
    (cl-assert (< (la-norm1-m (la-sub-mm (car (la-hermite [[1 1 2] [0 1 1] [1 2 3]])) [[1 0 1] [0 1 1] [0 0 0]])) 1e-5) "error in la-hermite (H not correct)")
    (cl-assert (< (la-norm1-m (la-sub-mm (cadr (la-hermite [[1 1 2] [0 1 1] [1 2 3]])) [[1 -1 0] [0 1 0] [-1 -1 1]])) 1e-5) "error in la-hermite (P not correct)")
    (cl-assert (equal (la-rank [[1 2 3] [3 6 9]]) 1) "error in la-rank")
    (cl-assert (equal (la-rank [[1 2 3] [3 6 7]]) 2) "error in la-rank")
    (cl-assert (< (la-norm1-v (la-sub-vv (la-project-vv [3.0 4.0 5.0] [3.0 4.0 -5.0]) [0 0 0]))) "error in la-project-vv")
    (cl-assert (< (la-norm1-v (la-sub-vv (la-project-vv [3.0 4.0 5.0] [1.0 0.0 0.0]) [1.0 0 0]))) "error in la-project-vv")
    (cl-assert (< (la-norm1-m (la-sub-mm (car (la-gram-schmidt-qr m6)) [[.6 .8] [.8 -.6]])) 1e-5) "error in la-gram-schmidt-qr (Q not correct)")
    (cl-assert (< (la-norm1-m (la-sub-mm (cadr (la-gram-schmidt-qr m6)) [[5 2] [0 1]])) 1e-5) "error in la-gram-schmidt-qr (R not correct)")
    (cl-assert (equal (la-eigenvalue [[2.0 6.0] [2.0 3.0]] 1e-50) '([6.0 -1.0] . [0.0 0.0])) "error in la-eigenvalue")
    (cl-assert (equal (la-eigenvalue [[0.0 -1.0] [1.0 0.0]]) '([0.0 0.0] . [1.0 -1.0])) "error in la-eigenvalue")
    (cl-assert (equal (la-eigenvalue [[1.0 1.0] [0.0 1.0]]) '([1.0 1.0] . [0.0 0.0])) "error in la-eigenvalue")
    
    )
  "test complete successfully"
  )

(message (la--test-function))

;;; test-linear-algebra ends here
