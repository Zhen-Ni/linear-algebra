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
	 )
	
    (cl-assert (equal (la-shape 9) nil) t "error in la-shape")
    (cl-assert (equal (la-shape v1) '(2)) t "error in la-shape")
    (cl-assert (equal (la-shape m2) '(3 2)) t "error in la-shape")
    (cl-assert (equal (la-transpose m2) m2t) t "error in la-transpose")
    (cl-assert (equal (la-identity 4) m3) t "error in la-identity")
    (cl-assert (equal (la-kron [[1 2] [3 4]] [[1 2] [3 4]]) [[1 2 2 4] [3 4 6 8] [3 6 4 8] [9 12 12 16]]) "error in la-kron")
    
    (cl-assert (equal (la-cwise-vv * v1 [4 5]) [4 10]) "error in la-cwise-vv")
    (cl-assert (equal (la-cwise-mm + m1 m4) [[6 8] [10 12]]) "error in la-cwise-mm")
    (cl-assert (equal (la-mul-vv v1 v2) (la-cwise-vv * v1 v2)) "error in la-mul-vv")
    (cl-assert (equal (la-add-mm m1 m4) (la-cwise-mm + m1 m4)) "error in la-add-mm")
    (cl-assert (equal (la-dot-vv [1 2 3] [3 2 1]) 10) "error in la-dot-vv")
    (cl-assert (equal (la-dot-mv m1 [2 1]) [4 10]) "error in la-dot-vm")
    (cl-assert (equal (la-dot-vm [2 1] (la-transpose m1)) [4 10]) "error in la-dot-mv")
    (cl-assert (equal (la-dot-mm m5 [[1] [2] [3]]) [[14] [32]]) "error in la-dot-mm")

    (cl-assert (equal (la-norm1-m [[1 2] [3 -4]]) 10) "error in la-norm1-v")
    (cl-assert (< (la-norm1-v (la-sub-vv (la-solve-mv m1 [4 10]) [2 1])) 1e-5) "error in la-solve-mv")
    (cl-assert (< (la-norm1-m (la-sub-mm (la-dot-mm m1 (la-inv m1)) (la-identity 2))) 1e-5) "error in la-solve-mv")
    (cl-assert (< (la-norm1-m (la-sub-mm (car (la-hermite [[1 1 2] [0 1 1] [1 2 3]])) [[1 0 1] [0 1 1] [0 0 0]])) 1e-5) "error in la-hermite (H not correct)")
    (cl-assert (< (la-norm1-m (la-sub-mm (cadr (la-hermite [[1 1 2] [0 1 1] [1 2 3]])) [[1 -1 0] [0 1 0] [-1 -1 1]])) 1e-5) "error in la-hermite (P not correct)")


    
    )
  "test complete"
  )

(message (la--test-function))

;;; test-linear-algebra ends here
