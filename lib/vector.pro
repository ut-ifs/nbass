Function scalar_pr_cyl,a,b
  ;calculates scalar product of two vectors in cylindrical coordinates
  ;transformation to cartezian [X,Y,Z], X along beam port R - right handed
  a_cart=[a(0)*cos(a(1)),a(0)*sin(a(1)),a(2)]
  b_cart=[b(0)*cos(b(1)),b(0)*sin(b(1)),b(2)]
  return,total(a_cart*b_cart)
end

Function scalar_pr_cart,a,b
  ;calculates scalar product of two vectors in cartesian coordinates
  return,total(a*b)
end

Function ang_2vect_cyl,a,b
  ;calculates angle between two vectors in cylindrical coordinates
  ang=acos(scalar_pr_cyl(a,b)/sqrt(scalar_pr_cyl(a,a)*scalar_pr_cyl(b,b)))
  return, ang
end

Function ang_2vect_cart,a,b
  ;calculates angle between two vectors in cartesian coordinates
  return,acos(scalar_pr_cart(a,b)/sqrt(scalar_pr_cart(a,a)*scalar_pr_cart(b,b)))
end

Function ang_2vect_norm,a,b
  ;calculates angle between two unit vectors in cartesian coordinates
  return,acos(scalar_pr_cart(a,b))
end

Function vect_from_2pts_cyl,A,B
  ;calculated scalar product of two vectros in cylindrical coordinates
  ;transformation to cartezian [X,Y,Z], X along beam port R - right handed
  A_cart=[A(0)*cos(A(1)),A(0)*sin(A(1)),A(2)];point A in cartesian
  B_cart=[B(0)*cos(B(1)),B(0)*sin(B(1)),B(2)];point B in cartesian
  AB_vect_cart=B_cart-A_cart
  AB_vect_cyl=[sqrt(AB_vect_cart(0)^2.0+AB_vect_cart(1)^2.0),atan(AB_vect_cart(1),AB_vect_cart(0)),AB_vect_cart(2)]
  return,AB_vect_cyl
end

Function cyl2cart,A
  return,[A(0)*cos(A(1)),A(0)*sin(A(1)),A(2)];point A in cartesian
end

Function cart2cyl,A
  return,[sqrt(A(0)^2.0+A(1)^2.0),atan(A(1),A(0)),A(2)]
end

Function norm_vect,A
  return,A/sqrt(total(a*a))
end

pro line_intersection,A1,A2,B1,B2,C1,C2,midpoint
  ;get the shortest line segment C between two skew lines A and B and find the midpoint
  R13 = A1-B1
  R43 = B2-B1
  R21 = A2-A1
  d1343 = scalar_pr_cart(R13,R43)
  d4321 = scalar_pr_cart(R43,R21)
  d1321 = scalar_pr_cart(R13,R21)
  d4343 = scalar_pr_cart(R43,R43)
  d2121 = scalar_pr_cart(R21,R21)
  denom = d2121*d4343 - d4321*d4321
  mua = (d1343*d4321 - d1321*d4343)/denom
  mub = (d1343 + d4321*mua)/d4343
  C1 = A1 + mua*(A2-A1)
  C2 = B1 + mub*(B2-B1)
  midpoint = (C1 + C2)/2
end
