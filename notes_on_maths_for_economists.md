--- 
title: "Notes on maths for economists"
author: "Juan Alvaro Díaz Raimond Kedilhac"
date: "2024-12-09"
site: bookdown::bookdown_site
documentclass: book
bibliography:
- book.bib
- packages.bib
description: |
  This is a minimal example of using the bookdown package to write a book.
  set in the _output.yml file.
  The HTML output format for this example is bookdown::gitbook,
link-citations: yes
github-repo: "rstudio/bookdown-demo"
---

# About  

Welcome to this **maths for economists study guide**, a resource designed to support undergraduate students in their journey to understanding maths. These personal notes are based on lectures attended at **Universidad Panamericana**, alongside additional materials such as lecture notes, websites, books, and other sources.  

While this guide is crafted to offer a quick and accessible approach to the subject, it is strongly recommended that students consult the original resources for a more comprehensive understanding.  

## Purpose  

This site serves as a developing repository of knowledge, aimed at supplementing students' studies. Please note that the material is continually updated, so your patience is appreciated as we work on corrections and add new content.  

## Feedback  

Your input is valuable! If you have comments, recommendations, or corrections, please reach out to us at **0251520@up.edu.mx**. We aim to address feedback and update the content promptly.  

## Support  

If you find this project helpful and wish to support it, you can make a donation at [buymeacoffee.com/jadrk040507](https://buymeacoffee.com/jadrk040507). Your generosity is greatly appreciated!  

Thank you for visiting, and happy learning!  




<!--chapter:end:index.Rmd-->

# Linear Algebra

## The Basic Vector Space \( \mathbb{R}^n \)

The quintessential vector space is \( \mathbb{R}^n \), which is defined as the set of all \( n \)-tuples of real numbers. As we will see in the next section, what makes this set a vector space is the ability to add two \( n \)-tuples to obtain another \( n \)-tuple:

$$
(v_1, \ldots, v_n) + (w_1, \ldots, w_n) = (v_1 + w_1, \ldots, v_n + w_n)
$$

and to multiply each \( n \)-tuple by a real number \( \lambda \):

$$
\lambda (v_1, \ldots, v_n) = (\lambda v_1, \ldots, \lambda v_n)
$$

In this way, each \( n \)-tuple is commonly referred to as a vector, and the real numbers \( \lambda \) are known as scalars. When \( n = 2 \) or \( n = 3 \), this reduces to vectors in the plane and in space, which most of us learned in high school.

The natural relationship from \( \mathbb{R}^n \) to \( \mathbb{R}^m \) is established through matrix multiplication. We write a vector \( \mathbb{x} \in \mathbb{R}^n \) as a column vector:

$$
\mathbb{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix}
$$

Similarly, we can write a vector in \( \mathbb{R}^m \) as a column vector with \( m \) entries. Let \( A \) be an \( m \times n \) matrix:

$$
A = \begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

Then, the product \( A \mathbb{x} \) is the \( m \)-tuple:

$$
A \mathbb{x} = \begin{pmatrix}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n \\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n \\
\vdots \\
a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n
\end{pmatrix}
$$

For any two vectors \( \mathbb{x} \) and \( \mathbb{y} \) in \( \mathbb{R}^n \) and any two scalars \( \lambda \) and \( \mu \), the following property holds:

$$
A (\lambda \mathbb{x} + \mu \mathbb{y}) = \lambda A \mathbb{x}  + \mu A \mathbb{y}
$$

In the next section, we will use the linearity of matrix multiplication to motivate the definition of a linear transformation between vector spaces. Now, let's relate all this to solving a system of linear equations. Suppose we are given numbers \( b_1, \ldots, b_m \) and numbers \( a_{ij}, \ldots, a_{mn} \). Our goal is to find \( n \) numbers \( x_1, \ldots, x_n \) that solve the following system of linear equations:

$$
\begin{aligned}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n &= b_1 \\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n &= b_2 \\
&\vdots \\
a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n &= b_m
\end{aligned}
$$

Calculations in linear algebra often reduce to solving a system of linear equations. When there are only a few equations, we can find the solutions manually, but as the number of equations increases, the calculations quickly become less about pleasant algebraic manipulations and more about keeping track of many small individual details. In other words, it is an organizational problem.

We can write:

$$
\mathbb{b} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}
$$

and our unknowns as:

$$
\mathbb{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix}
$$

Then, we can rewrite our system of linear equations in the more visually appealing form:

$$
A \mathbb{x} = \mathbb{b}
$$

When \( m > n \) (when there are more equations than unknowns), we generally do not expect solutions. For example, when \( m = 3 \) and \( n = 2 \), this corresponds geometrically to the fact that three lines in a plane generally do not have a common intersection point. When \( m < n \) (when there are more unknowns than equations), we generally expect there to be many solutions. In the case where \( m = 2 \) and \( n = 3 \), this corresponds geometrically to the fact that two planes in space generally intersect in an entire line. Much of the machinery of linear algebra deals with the remaining case when \( m = n \).

Therefore, we want to find the \( n \times 1 \) column vector \( \mathbb{x} \) that solves \( A \mathbb{x} = \mathbb{b} \), where \( A \) is a given \( n \times n \) matrix and \( \mathbb{b} \) is a given \( n \times 1 \) column vector.

Suppose the square matrix \( A \) has an inverse matrix \( A^{-1} \) (which means that \( A^{-1} \) is also \( n \times n \) and, more importantly, that \( A A^{-1} = I \), where \( I \) is the identity matrix). Then our solution will be:

$$
\mathbb{x} = A^{-1} \mathbb{b}
$$

because:

$$
A \mathbb{x} = A (A^{-1} \mathbb{b}) = I \mathbb{b} = \mathbb{b}
$$

Thus, solving our system of linear equations reduces to understanding when the \( n \times n \) matrix \( A \) has an inverse. (If an inverse matrix exists, then algorithms exist for its calculation). The key theorem of linear algebra, which is stated in section six, is essentially a list of many equivalencies for when an \( n \times n \) matrix has an inverse, and it is therefore crucial to understanding when a system of linear equations can be solved.




## Vector Spaces and Linear Transformations

The abstract approach to studying systems of linear equations begins with the notion of a vector space.

::: {.definition}  
A set \( V \) is a vector space over the real numbers \( \mathbb{R} \) if there are two operations:

1. Scalar multiplication: For every \( a \in \mathbb{R} \) and \( \mathbb{v} \in V \), there is an element \( a \cdot \mathbb{v} \in V \), denoted by \( a \mathbb{v} \), which satisfies the properties of scalar multiplication.

2. Vector addition: For every \( \mathbb{v}, \mathbb{w} \in V \), there is an element \( \mathbb{v} + \mathbb{w} \in V \), denoted by \( \mathbb{v} + \mathbb{w} \), which satisfies the properties of vector addition.

These operations must satisfy the following properties:

a) **Additive identity**: There exists an element \( \mathbb{0} \in V \) (the zero vector) such that \( \mathbb{0} + \mathbb{v} = \mathbb{v} \) for all \( \mathbb{v} \in V \).

b) **Additive inverse**: For each \( \mathbb{v} \in V \), there exists an element \( -\mathbb{v} \in V \) such that \( \mathbb{v} + (-\mathbb{v}) = \mathbb{0} \).

c) **Commutativity**: For all \( \mathbb{v}, \mathbb{w} \in V \), it holds that \( \mathbb{v} + \mathbb{w} = \mathbb{w} + \mathbb{v} \).

d) **Distributivity of scalar multiplication over vector addition**: For all \( a \in \mathbb{R} \) and for all \( \mathbb{v}, \mathbb{w} \in V \), we have \( a(\mathbb{v} + \mathbb{w}) = a\mathbb{v} + a\mathbb{w} \).

e) **Distributivity of scalar multiplication over scalar addition**: For all \( a, b \in \mathbb{R} \) and for all \( \mathbb{v} \in V \), it holds that \( (a + b)\mathbb{v} = a\mathbb{v} + b\mathbb{v} \).

f) **Compatibility of scalar multiplication with field multiplication**: For all \( a, b \in \mathbb{R} \) and all \( \mathbb{v} \in V \), it holds that \( a(b\mathbb{v}) = (ab)\mathbb{v} \).

g) **Multiplicative identity**: For all \( \mathbb{v} \in V \), it holds that \( 1 \cdot \mathbb{v} = \mathbb{v} \), where 1 is the multiplicative identity in \( \mathbb{R} \).
:::

Real numbers can be replaced by complex numbers and, indeed, by any field.

As a matter of notation, and to be consistent with common usage, the elements of a vector space are called vectors, and the elements of \( \mathbb{R} \) (or any field being used) are called scalars. It is worth noting that the space \( \mathbb{R}^n \) mentioned in the previous section satisfies these conditions.

The natural map between vector spaces is that of a linear transformation.

::: {.definition} A linear transformation \( T : V \to W \) is a function from a vector space \( V \) to a vector space \( W \) such that for any real numbers \( a_1 \) and \( a_2 \), and for any vectors \( v_1 \) and \( v_2 \) in \( V \), it holds that

$$
T(a_1 v_1 + a_2 v_2) = a_1 T(v_1) + a_2 T(v_2).
$$

Matrix multiplication from \( \mathbb{R}^n \) to \( \mathbb{R}^m \) provides an example of a linear transformation.
:::

::: {.definition} A subset \( U \) of a vector space \( V \) is a subspace of \( V \) if \( U \) itself is a vector space.

In practice, it is usually easy to determine if a subset of a vector space is, in fact, a subspace, using the following proposition, whose proof is left to the reader.
:::

::: {.proposition} A subset \( U \) of a vector space \( V \) is a subspace of \( V \) if \( U \) is closed under addition and scalar multiplication.

Given a linear transformation \( T : V \to W \), there are naturally occurring subspaces both in \( V \) and \( W \).
:::

::: {.definition} If \( T : V \to W \) is a linear transformation, then the kernel of \( T \) is:

$$
\ker(T) = \{ v \in V : T(v) = 0 \}
$$

and the image of \( T \) is

$$
\text{Im}(T) = \{ w \in W : \text{there exists } v \in V \text{ such that } T(v) = w \}.
$$

The kernel is a subspace of \( V \), since if \( v_1 \) and \( v_2 \) are two vectors in the kernel and \( a \) and \( b \) are any two real numbers, then

$$
T(a v_1 + b v_2) = a T(v_1) + b T(v_2) = a \cdot 0 + b \cdot 0 = 0.
$$
:::

Similarly, we can prove that the image of \( T \) is a subspace of \( W \).

If the only vector spaces that existed were column vectors in \( \mathbb{R} \), then even this level of abstraction would be trivial. However, this is not the case.

Here we consider just one example. Let \( C^*[0,1] \) be the set of all real-valued functions with domain on the unit interval \( [0,1] \):

$$
f : [0,1] \to \mathbb{R}
$$

such that the \( k \)-th derivative of \( f \) exists and is continuous. Since the sum of any two such functions and a multiple of any of these functions by a scalar will remain in \( C^*[0,1] \), we have a vector space. Although we will define the dimension officially in the next section, \( C^*[0,1] \) will have infinite dimension (and thus will definitely not be some \( \mathbb{R}^n \)). We can see the derivative as a linear transformation from \( C^*[0,1] \) to those functions with one derivative less, \( C^{k-1}[0,1] \):

$$
\frac{d}{dx} : C^*[0,1] \to C^{k-1}[0,1].
$$

The kernel of this transformation consists of those functions whose \( k \)-th derivative is zero, that is, constant functions.

Now consider the differential equation

$$
f'' + 3f' + 2f = 0.
$$

Let \( T \) be the linear transformation:

$$
T : C^*[0,1] \to C^*[0,1]
$$

defined by

$$
T(f) = f'' + 3f' + 2f.
$$

The problem of finding a solution \( f(x) \) to the original differential equation can now be translated into finding an element in the kernel of \( T \). This suggests the possibility (which is indeed true) that the language of linear algebra can be used to understand solutions to (linear) differential equations.




## Bases, Dimension, and Linear Transformations as Matrices

Our next goal is to define the dimension of a vector space.

::: {.definition}
A set of vectors \((v_1, \dots, v_n)\) forms a basis for the vector space \(V\) if, for any vector \(v \in V\), there exist unique scalars \(a_1, \dots, a_n \in \mathbb{R}\) such that 

\[
v = a_1v_1 + \dots + a_nv_n.
\]
:::

::: {.definition}
The dimension of a vector space \(V\), denoted as \(\text{dim}(V)\), is the number of elements in a basis.
:::

It is not obvious that the number of elements in a basis will always be the same, regardless of the chosen basis. To ensure that the definition of the dimension of a vector space is well-defined, we need the following theorem (which we will not prove):

::: {.definition}
All bases of a vector space \(V\) have the same number of elements.
:::

For \(\mathbb{R}^n\), the usual basis is \(\{(1, 0, \dots, 0), (0, 1, 0, \dots, 0), \dots, (0, \dots, 0, 1)\}\). Therefore, \(\mathbb{R}^n\) has dimension \(n\). If this were not the case, the previous definition of dimension would be incorrect and we would need another. This is an example of the principle mentioned in the introduction. We have an intuitive understanding of what the dimension should mean for specific examples: a line should be one-dimensional, a plane two-dimensional, and three-dimensional space. We then formulate a precise definition. If this definition gives the "correct answer" for our three already understood examples, we are somewhat confident that the definition has captured what dimension means in this case. We can then apply the definition to examples where our intuitions fail.

Linked to the idea of a basis is:

::: {.definition}
The vectors \((v_1, \dots, v_n)\) in a vector space \(V\) are linearly independent if, whenever \(a_1v_1 + \dots + a_nv_n = 0\), it must be the case that the scalars \(a_1, \dots, a_n\) are all zero. Intuitively, a set of vectors is linearly independent if they all point in different directions. A basis, therefore, consists of a set of linearly independent vectors that span the vector space, where "span" means:
:::

::: {.definition}
A set of vectors \((v_1, \dots, v_n)\) spans the vector space \(V\) if, for any vector \(v \in V\), there exist scalars \(a_1, \dots, a_n \in \mathbb{R}\) such that

\[
v = a_1v_1 + \dots + a_nv_n.
\]
:::

Our next goal is to show how all linear transformations \(T: V \to W\) between finite-dimensional spaces can be represented as matrix multiplication, provided that we fix bases for the vector spaces \(V\) and \(W\).

First, we fix a basis \(\{v_1, \dots, v_n\}\) for \(V\) and a basis \(\{w_1, \dots, w_m\}\) for \(W\). Before examining the linear transformation \(T\), we need to show how each element of the \(n\)-dimensional space \(V\) can be represented as a column vector in \(\mathbb{R}^n\) and how each element of the \(m\)-dimensional space \(W\) can be represented as a column vector in \(\mathbb{R}^m\). Given any vector \(v \in V\), by the definition of a basis, there exist unique real numbers \(a_1, \dots, a_n\) such that:

\[
v = a_1v_1 + \dots + a_nv_n.
\]

Thus, we represent the vector \(v\) with the column vector:

\[
\begin{pmatrix} a_1 \\ a_2 \\ \vdots \\ a_n \end{pmatrix}.
\]

Similarly, for any vector \(w \in W\), there exist unique real numbers \(b_1, \dots, b_m\) such that:

\[
w = b_1w_1 + \dots + b_mw_m.
\]

Here, we represent \(w\) as the column vector:

\[
\begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}.
\]

It is important to note that we have established a correspondence between the vectors in \(V\) and \(W\) and the column vectors in \(\mathbb{R}^n\) and \(\mathbb{R}^m\), respectively. More technically, we can prove that \(V\) is isomorphic to \(\mathbb{R}^n\) (which means there is a one-to-one and onto linear transformation from \(V\) to \(\mathbb{R}^n\)) and that \(W\) is isomorphic to \(\mathbb{R}^m\), although it must be emphasized that the real correspondence only exists after a basis has been chosen (which means that, while the isomorphism exists, it is not canonical; this is an important aspect because, in practice, we are often not provided with a basis).

Now, we want to represent a linear transformation \(T: V \to W\) as a matrix \(A\) of size \(m \times n\). For each basis vector \(v_i\) in the vector space \(V\), \(T(v_i)\) will be a vector in \(W\). Therefore, there will be real numbers \(a_{ij}, \dots, a_{im}\) such that:

\[
T(v_i) = a_{ij}w_1 + \dots + a_{im}w_m.
\]

We want to see that the linear transformation \(T\) corresponds to the matrix \(A\):

\[
A = \begin{pmatrix} a_{11} & a_{21} & \dots & a_{m1} \\ a_{12} & a_{22} & \dots & a_{m2} \\ \vdots & \vdots & \ddots & \vdots \\ a_{1n} & a_{2n} & \dots & a_{mn} \end{pmatrix}.
\]

Given any vector \(v \in V\), with \(v = a_1v_1 + \dots + a_nv_n\), we have:

\[
T(v) = a_1T(v_1) + \dots + a_nT(v_n) = a_1(a_{11}w_1 + \dots + a_{1m}w_m) + \dots + a_n(a_{n1}w_1 + \dots + a_{nm}w_m).
\]

Under the correspondences of the vector spaces with the respective column spaces, this can be seen as the matrix multiplication of \(A\) by the column vector corresponding to the vector \(v\):

\[
A \begin{pmatrix} a_1 \\ a_2 \\ \vdots \\ a_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}.
\]

It is important to note that if \(T: V \to V\) is a linear transformation of a vector space onto itself, then the corresponding matrix will be \(n \times n\), i.e., a square matrix.

Given different bases for the vector spaces \(V\) and \(W\), the matrix associated with the linear transformation \(T\) will change. A natural problem is to determine when two matrices actually represent the same linear transformation, but under different bases. This will be the subject of section seven.



<!--chapter:end:07-notes_econometrics.Rmd-->



<!--chapter:end:07-references.Rmd-->

