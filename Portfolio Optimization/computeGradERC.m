function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ;  

  if(size(x,1)==1)
     x = x';
  end
  
  % TO CALCULATE THE GRADIENT OF OBJECTIVE FUNCTION
  
  % calculate number of pairs (RC(i) - RC(j))
  num_pairs = (n*n - n)/2;
  
  % create a matrix to store derivatives of each pair in respect to each of
  % the w(i) 
  der_pair = zeros(n, num_pairs);
  
  row = Q*x;
  RC = x .* (Q*x) ;
  pair_count = 0;
  
  % do the loop to get the derivatives of each pair in respect to each of
  % w(i) 
  for i = 1:n;
     for j = i+1:n;
%         for h = 1:num_pairs;
            pair_i_j = RC(i) - RC(j);
            pair_count = pair_count + 1;
            for k = 1:n;
                if k == i;
                    % multiply by 2 because the occurence of one pair [RC(i) - RC(j)] is 2
                    der_pair(k,pair_count) = 2*2*pair_i_j*(row(i) + x(i)*Q(i,i) - x(j)*Q(i,j));
                
                elseif k == j;
                    der_pair(k,pair_count) = 2*2*pair_i_j*(x(i)*Q(i,j) - row(j) - x(j)*Q(i,j));
                
                else 
                    der_pair(k,pair_count) = 2*2*pair_i_j*(x(i)*Q(i,k) - x(j)*Q(j,k));
                end
            end
     end
  end
  
  % the sum of each of the row is the derivatives in respect to w(i) of the
  % function 
  
  gval = sum(der_pair,2); 
 
end
