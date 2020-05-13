function result=trianglenumbers(n)
% Generates the triangle number, which is the factorial, but addition.
% That is, the sum of natural numbers
% Returns the total number of elements in the triangle
    result=[];
    sum=0;
    for i = 1:n
        sum=sum+i;
        result=[sum];
    end
end