function TitleScreen(w, center, message)

% start screen for fmri experiments.  may add flag to say press button or
% that it will start with the scanner.  takes 3 arguments
% w : screen pointer
% center: position on the screen to center the text around
% message: text to display

 Screen('TextSize',w,40);
%  get size of text and display on center of screen
width = RectWidth(Screen('TextBounds',w,message))
 Screen('DrawText',w,message,center(1)-width/2, center(2),[0 0 0]);
 Screen('Flip',w);

end