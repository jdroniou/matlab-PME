% Transform between the real function u, and the vector of unknowns X used in computations
	%		m>1: X is the solution 'u' at the vertices
  %   m<1: X is |u|^(m-1)u at the vertices

function a=transform(b,direction)

global m;
a=zeros(size(b));
sb = sign(b);

if (direction=='to_unknowns')
  % b is the real function, a are the unknowns
  if (m>1)
    a = b;
  else
    a = abs(b).^m.*sb;  
  end

else
  % b is the unknown, a is the function
  if (m>1)
    a = b;
  else
    a = abs(b).^(1/m).*sb;
  end

end


