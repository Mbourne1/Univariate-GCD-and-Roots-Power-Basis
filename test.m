
a = ['red' 'blue']
b = ['car' 'bike']
[A,B] = meshgrid(a,b);

c=cat(2,A',B');

d=reshape(c,[],2)

