function [linkWeights] = dijkstra(A, startNode, maxDest) 
%
% Modified: 03/03/2003 by Deepayan Chakrabarti
% Obtained from: http://www.owlnet.rice.edu/~blakeb/vlsi/overallalgo.html
%
% Params: A is the matrix (currently we require it to be square)
%         startNode -> the "source" node
%         maxDest -> Find shortest paths to how many destinations from the
%                    source. This lets the algo terminate
%                    w/o going through the entire connected component.
%                    Set to -1 to do the entire component.
%
% Moffett Stephen, Jacob Rios, Angel Sun and Blake Borgeson 
% 9/28/01 
% ELEC 422 Group D 
% 
% Solving the Shortest Path Problem using Dijkstra's Algorithm 
% - This function takes in a square matrix as input. The matrix 
% represents the weights of edges in a network of n vertices. A zero 
% indicates that vertices are not adjacent. The function returns 
% an n-1 x 3 matrix, sp, that shows the shortest total path length 
% to each vertex from the starting vertex (in the middle column), 
% and the preceding vertex in the minimum path (in the right column). 
% 
% DIJKSTRA find the shortest paths from the starting vertex to 
% each other vertex in the network. 
% 
% INPUT: 
% A = input matrix 
% 
% OUTPUT: 
% sp = shortest paths to each vertex from starting vertex 

[m,n] = size(A); 
y = [zeros(1,n); A(startNode,1:n); ones(1,n)*startNode]; % shortest path 
y(1,startNode) = 1;
addedNodes = [zeros(1,n-1)]; % contains the nodes in the order in which they
                             % get added.
			     
linkWeights = sparse(m,n);  % contains the link weights for each link
			     
% status (volatile memory) 

p = find(y(1,1:n)==0); % p = vector of non-terminated vertices 
iter = 1; 

while ~isempty(p) 
% fprintf('\n iter #%d: y = \n ',iter); disp(y); 
% pause; 
x = find(y(2,p)>0); % x = vector of indices referring to nonzero entries 
if isempty(x), break, end

% Step 2 
[a,K] = min(y(2,p(x))); % a = smallest nonzero entry 

% K = index of x referring to a 
o = p(x); 
J = o(K);     % J = index of A referring to a 
y(1,J) = 1;   % changes termination bit to 1 
p = p(find(p ~= J)); % update p 
addedNodes(iter) = J;

if ~isempty(p) 
z = find(A(J,p)>0); % z = vector of vertices adjacent to J 
r = p(z); % r = indices (wrt A) of non-terminated adjacent vertices 
w = A(J,r) + y(2,J); % distance to all vertices via J 
temp1 = y(2,r); 
temp2 = y(3,r); 
helper1 = y(2,r) > w; % if w is less than previous 

% shortest path ... 
y(3,r) = y(3,r) + (J - y(3,r)).*helper1; % update y 
y(2,r) = min(y(2,r), w); 
helper2 = temp1 == 0; % if shortest path was zero ... 
y(2,r) = y(2,r) + (w - y(2,r)).*helper2; % update y 
y(3,r) = y(3,r) + (J - y(3,r)).*helper2; 
% y stays unchanged if w > y(2,r-1) 

end 
iter = iter + 1; 
if ((maxDest>=1)  & (iter>maxDest)),
	break;
end

end 
%sp = [(2:n)',y(2,:)',y(3,:)'] 

%binform = [dec2bin(sp(:,1)) dec2bin(sp(:,2)) dec2bin(sp(:,3))]; 

% Set linkWeights. Start from the last addedNodes, and seek backwards.
for count = iter-1 : -1 : 1,
	tmp = addedNodes(count);
%	fprintf('count = %d, last node = %d, tmp = %d\n',count,y(3,tmp),tmp);
	linkWeights(tmp,y(3,tmp)) = linkWeights(tmp,y(3,tmp)) + sum(linkWeights(:,tmp)) + 1;
end

