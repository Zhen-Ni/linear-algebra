;;; linear-algebra.el --- linear algebra
;;; Commentary:
;;; This package provides the basics operations for linear algebra.



;;; Code:

(require 'cl-macs)
(require 'seq)


;;; Basic storage operations.


(defun la-vector (n &optional init)
  "Return a vector of size N, with values initialized to INIT."
  (make-vector n (or init 0))
  )

(defun la-matrix (n1 n2 &optional init)
  "Return a matrix of size N1 × N2, with values initialized to INIT."
  (let ((res (make-vector n1 nil)))
    (dotimes (i n1 res)
      (setf (aref res i) (la-vector n2 init))
      )
    )
  )

(defun la-identity (n)
  "Return an identity matrix of size N × N."
  (let ((res (la-matrix n n 0)))
    (dotimes (i n res)
      (setf (aref (aref res i) i) 1)
      )
    )
  )

(defun la-row (m n)
  "Return a vector which the elements are the N th row of matrix M."
  (aref m n)
)
  
(defun la-col (m n)
  "Return a vector which the elements are the N th column of matrix M."
  (let* ((size (car (la-shape m)))
	 (res (make-vector size nil))
	 )
    (dotimes (i size res)
      (setf (aref res i) (aref (aref m i) n))
      )
    )
  )

(defun la-shape (x)
  "Get the shape of matrix or vector or scalar X."
  (if (arrayp x)
      (let ((n1 (length x)))
	(if (equal (length x) 0)
	    '(0)
	  (let ((e (aref x 0)))
	    (if (arrayp e)
		`(,(length x) ,(length e))
	      `(,(length x))
	      )
	    )
	  )
	)
    ;; else branch of the outmost if
    `()
    )
  )


(defun la-transpose (m)
  "Transpose matrix M."
  (let* ((n1 (length m))
	 (n2 (length (aref m 0)))
	 (res (make-vector n2 nil))
	 )
    (dotimes (i n2 res)
      (let ((row (setf (aref res i) (make-vector n1 nil))))
	(dotimes (j n1)
	  (setf (aref row j) (aref (aref m j) i))
	  )
	)
      )
    )
  )


(defun la-kron (a b)
  "Kronecker product̨̨ of matrix A and B."
  (let* ((shape1 (la-shape a))
	 (shape2 (la-shape b))
	 (n1 (car shape1))
	 (n2 (cadr shape1))
	 (n3 (car shape2))
	 (n4 (cadr shape2))
	 (nx (* n1 n3))
	 (ny (* n2 n4))
	 (res (la-matrix nx ny))
	 )
    (dotimes (i nx res)
      (dotimes (j ny)
	(setf (aref (aref res i) j)
	      (* (aref (aref a (/ i n3)) (/ j n4))
		 (aref (aref b (% i n3)) (% j n4))
		 ))))
    )
  )

;;; Arithmetical operations.

;; Macro to avoid repeated code.
(defmacro la-cwise-vv (op v1 v2)
  "Do element-wise operation OP between vectors V1 and V2."
  `(let* ((v1 ,v1)
	  (v2 ,v2)
	  (n (length v1))
	 (res (make-vector n nil)))
    (cl-assert (equal n (length v2)) t "shape not match")
    (dotimes (i n res)
      (setf (aref res i) (,op (aref v1 i) (aref v2 i)))
      )
    )
  )

(defun la-add-vv (v1 v2)
  "Element-wise add of vectors V1 and V2."
  (la-cwise-vv + v1 v2)
  )

(defun la-sub-vv (v1 v2)
  "Element-wise substraction of vectors V1 and V2."
  (la-cwise-vv - v1 v2)
  )

(defun la-mul-vv (v1 v2)
  "Element-wise multipilication of vectors V1 and V2."
  (la-cwise-vv * v1 v2)
  )

(defun la-div-vv (v1 v2)
  "Element-wise division of vectors V1 and V2."
  (la-cwise-vv / v1 v2)
  )

;; Macro to avoid repeated code.
(defmacro la-cwise-mm (op m1 m2)
  "Do element-wise operation OP between matrixes M1 and M2."
  `(la-cwise-vv (lambda (v1 v2) (la-cwise-vv ,op v1 v2)) ,m1 ,m2)
  )

(defun la-add-mm (m1 m2)
  "Element-wise add of matrixes M1 and M2."
  (la-cwise-mm + m1 m2)
  )

(defun la-sub-mm (m1 m2)
  "Element-wise substraction of matrixes M1 and M2."
  (la-cwise-mm - m1 m2)
  )

(defun la-mul-mm (m1 m2)
  "Element-wise multipilication of matrixes M1 and M2."
  (la-cwise-mm * m1 m2)
  )

(defun la-div-mm (m1 m2)
  "Element-wise division of matrixes M1 and M2."
  (la-cwise-mm / m1 m2)
  )

(defun la-dot-vv (v1 v2)
  "Inner product of vectors V1 and V2."
  (seq-reduce '+ (la-cwise-vv * v1 v2) 0)
  )


(defun la-dot-mv (m v)
  "Matrix product of matrix M and vector V."
  (let* ((n (length m))
	 (res (make-vector n nil)))
    (dotimes (i n res)
      (setf (aref res i) (la-dot-vv (aref m i) v))
      )
    )
  )

(defun la-dot-vm (v m)
  "Matrix product of vector V and Matrix M."
  (la-dot-mv (la-transpose m) v)
  )

(defun la-dot-mm (m1 m2)
  "Matrix product of matrixes M1 and M2."
  (let* ((m2t (la-transpose m2))
	 (n1 (length m1))
	 (n2 (length m2t))
	 (res (make-vector n2 nil)))
    (la-transpose
     (dotimes (i n2 res)
       (setf (aref res i) (la-dot-mv m1 (aref m2t i)))
       )
     )
    )
  )


;;; Linear algebra operations.

(defun la-norm1-v (v)
  "The L1-norm of vector V."
  (let ((res 0))
    (seq-doseq (x v)
      (setf res (+ res (abs x)))
      )
    res)
  )

(defun la-norm1-m (m)
  "The L1-norm of matrix M."
  (let ((res 0))
    (seq-doseq (v m)
      (setf res (+ res (la-norm1-v v)))
      )
    res)
  )

(defun la--exchange-row (m l1 l2)
  "Exchange the rows L1 and L2 of matrix M.
It is an inplace operation, the return value is non-sense."
  (let ((tmp (aref m l1)))
    (setf (aref m l1) (aref m l2))
    (setf (aref m l2) tmp)
    )
  )

(defun la--multiply-row (m l c)
  "Multiply row L of matrix M by scalar c.
It is an inplace operation, the return value is non-sense."
  (setf (aref m l)
	(cl-map 'vector (lambda (i) (* c i)) (aref m l)))
  )
  
(defun la--stack-row (m l1 l2 c)
  "Add row L1 of matrix M by the multipication of row L2 and scalar c.
It is an inplace operation, the return value is non-sense."
  (setf (aref m l1)
	(la-add-vv
	 (aref m l1)
	 (cl-map 'vector (lambda (i) (* c i)) (aref m l2))))
  )


(defun la--argmaxabs (vec)
  "Find the index of the max abs value in VEC."
  (let ((maxval 0)
	(idx))
    (dotimes (i (length vec))
      (if (< (abs maxval) (abs (elt vec i)))
	  (setf idx i maxval (elt vec i))
	)
      )
    `(,idx ,maxval)
    )
  )

(defun la--slice (vec start end)
  "Get a subset of vector VEC from element START to END."
  (let ((res (make-vector (- end start) nil)))
    (dotimes (i (- end start) res)
	     (setf (aref res i) (aref vec (+ start i)))
	     ))
  )

(defun la-solve-mm (A B)
  "Solve matrix function A · X = B."
  (let* ((m (copy-tree A t))
	 (res (copy-tree B t))
	 (shape (la-shape m))
	 (shape2 (la-shape res))
	 (n (car shape)))
    (cl-assert (equal n (cadr shape)) t "A should be square")
    (cl-assert (equal n (car shape2)) t "dimension not match")
    (dotimes (i n res)
      (let* ((row (la--slice (la-col m i) i n))
	     (idx-maxval (la--argmaxabs row)) ; find the maximum
	     (c))
	(cl-assert (car idx-maxval) t "singular matrix")
	;; Exchange max value to the top to reduce error.
	(setf c (+ i (car idx-maxval)))
	(la--exchange-row m c i) (la--exchange-row res c i)
	;; Make this line begin with 1.
	(setf row (la--slice (la-col m i) i n)
	      c (/ 1.0 (aref row 0)))
	(la--multiply-row m i c) (la--multiply-row res i c)
	;; Elimate the first element in the remaining rows.
	(setf row (la-col m i))
	(dotimes (j n)
	  (if (equal i j) ()
	    (setf c (- (aref row j)))
	    (la--stack-row m j i c) (la--stack-row res j i c)
	    ))
	)))
  )
  
(defun la-solve-mv (A b)
  "Solve matrix function A · X = B."
  (la-col (la-solve-mm A (la-transpose (vector b))) 0)
  )

(defun la-inv (m)
  "Inverse matrix M."
  (la-solve-mm m (la-identity (car (la-shape m))))
  )

(defun la-hermite (a &optional p)
  "Hermite decomposition of matrix A.
This function returns two matrixes, P and H, which satisfies
P · A = H, where H is the hermite form and P represents the
transformation.  The optional argument p is used for record the
transformation, which defaults to identity matrix."
  (let* ((h (copy-tree a t))
	 (shape (la-shape h))
	 (n1 (car shape))
	 (n2 (cadr shape))
	 (p (or (copy-tree p t) (la-identity n1)))
	 (n (seq-min shape))
	 )
    (dotimes (i n (list h p))
      ;; Find the maximum value in the row.
      (let* ((row (la--slice (la-col h i) i n1))
	     (idx-maxval (la--argmaxabs row))
	     (c))
	;; Do operations only when the row contains at least one value
	;; that is non-zero.
	(if (car idx-maxval)
	    (progn
	      ;; Exchange max value to the top to reduce error.
	      (setf c (+ i (car idx-maxval)))
	      (la--exchange-row h c i) (la--exchange-row p c i)
	      ;; Make this line begin with 1.
	      (setf row (la--slice (la-col h i) i n1)
		    c (/ 1.0 (aref row 0)))
	      (la--multiply-row h i c) (la--multiply-row p i c)
	      ;; Elimate the first element in the remaining rows.
	      (setf row (la-col h i))
	      (dotimes (j n)
		(if (equal i j) ()
		  (setf c (- (aref row j)))
		  (la--stack-row h j i c) (la--stack-row p j i c)
		  )))))))
  )



(provide 'linear-algebra)

;;; linear-algebra ends here
