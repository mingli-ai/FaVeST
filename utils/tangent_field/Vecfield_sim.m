function f = flow_rand(X)

rng(111);
N = size(X,2);
curl_free_dimx = rand(1,N);
curl_free_dimy = rand(1,N);
curl_free_dimz = rand(1,N);
T1 = [curl_free_dimx;curl_free_dimy;curl_free_dimz];
for i=1:3
    T(i,:)= T1(i,:)-sum(X.*T1).*T1(i,:);
end

curl_free = T';

diver_free_dimx = rand(1,N);
diver_free_dimy = rand(1,N);
diver_free_dimz = rand(1,N);
T2 = [diver_free_dimx;diver_free_dimy;diver_free_dimz];
for i=1:3
    T(i,:)= T2(i,:)-sum(X.*T2).*T2(i,:);
end

diver_free = T';

f = curl_free+diver_free;


end

