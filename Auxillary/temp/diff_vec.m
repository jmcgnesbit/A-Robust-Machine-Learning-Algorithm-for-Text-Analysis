function p = diff_vec(u,v)
% p = diff_vec(u,v)
% Returns a difference ratio of v, w.r.t. u.

%[M N]=size(u);
%p=(ones(1,M)*abs(u-v)*ones(N,1))/(ones(1,M)*abs(u)*ones(N,1));
% Original
p = norm(u - v) / norm(u);
