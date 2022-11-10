function f = landing_gear(q,s,r,r2,L1,L2,d1,d2,d3,d4,d5,d6)

f(1) = -(r+r2+s)*cos(q(1))+d4*sin(q(2))-d6;
f(2) = -(r+r2+s)*sin(q(1))-d4*cos(q(2))+d1;
f(3) = -d2+d5*sin(q(3))+L2*cos(q(4))+L1*sin(q(2));
f(4) = -d3+d5*cos(q(3))+L2*sin(q(4))-L1*cos(q(2))+d1;

end

