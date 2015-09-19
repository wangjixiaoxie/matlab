function key=KeyPressTestFunction
% Attaches KeyPress test function to the figure
%  T. C. O'Haver (toh@umd.edu),  Version 2, March 2008.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')

function ReadKey(obj,eventdata)
% When a key is pressed, interprets the key and calls
%  corresponding function.
key=get(gcf,'CurrentCharacter');

if ischar(key),
  switch double(key),
    case  28
        % Pans one point down when left arrow pressed.
        ipfnudgedown
    case 29
        % Pans one point up when right arrow pressed.
        ipfnudgeup
    case 30
        % Zooms one point up when up arrow pressed.
        ipfzoomup
    case 31
        % Zooms one point down when down arrow pressed.
        ipfzoomdown
    case 102
        % When 'f' key is pressed,
        % Fits and plots, same as clicking the 'Re-fit' slider
        ipfrefit(.1,gca)
    case 99
        % When 'c' key is pressed, user clicks graph to enter start positons, 
        % then fit is computed and graph re-drawn.
        % Same as clicking the 'Custom' slider
        ipfcustom    
    case 98
        % When 'b' key is pressed, user clicks graph to enter background points, 
        % then fit is computed and graph re-drawn.
        % Same as clicking the 'BG' slider
        ipfbackground
    case {49,50,51,52,53,54}
        % When a number key is pressed, sets the number of peaks to that
        % number.  Same as selecting with the % Peaks slider
        ipfpeaks(key-48,gca)
    case 103
        % When 'g' key is pressed, peak shape is set to
        % Gaussian. Same as using the Shape slider.
        ipfShape(1,gca)
     case 108
        % When 'l' key is pressed, peak shape is set to
        % Lorentzian. Same as using the Shape slider.
        ipfShape(2,gca)
     case 111
        % When 'o' key is pressed, peak shape is set to
        % Logistic. Same as using the Shape slider.
        ipfShape(3,gca)
     case 112
        % When 'p' key is pressed, peak shape is set to
        % Pearson. Same as using the Shape slider.
        ipfShape(4,gca)
     case 101
        % When 'e' key is pressed, peak shape is set to Exponentally-
        % broadened gaussian. Same as using the Shape slider.
        ipfShape(5,gca)

  end
end
    