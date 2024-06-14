function fig(s)
% FIG  Creates a new figure window with title.
%
% s: Title (optional string or number)

% James Taylor
% 06/08/1999

if nargin<1 | isempty(s)
  s=' ';
end

if ~ischar(s)
  s=num2str(s);
end

figure; zoom on
set(gcf, 'name', s)
set(gcf, 'NumberTitle', 'off')

% end of m-file