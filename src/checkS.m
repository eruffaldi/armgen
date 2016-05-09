

syms a d alpha theta dvar thetavar rx ry rz tx ty tz real;

Asr = rotz(theta+thetavar)*rotz(rz)*roty(ry)*rotx(rx)*trans([tx ty tz]);
psr = simplify(Asr(1:3,1:3)'*Asr(1:3,4))

Asp = trans([0 0 d+dvar])*rotz(rz)*roty(ry)*rotx(rx)*trans([tx ty tz]);
psp = simplify(Asp(1:3,1:3)'*Asp(1:3,4))

Aj = trans([0 0 d+dvar])*rotz(theta+thetavar)*trans([a 0 0])*rotx(alpha);
pj = simplify(Aj(1:3,1:3)'*Aj(1:3,4))