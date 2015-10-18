function[response]=subjectCResponse(w,rect,todisplay,drect,textDown);
%function that displays and records the subject's response as they type it;
%as the screen must be redrawn with each key press, the image that needs to
%be displayed while this happens must be passed in
%
%w is the window pointer
%rect is the size of the screen rect  (should be eliminated, don't need to
%pass this)
%todisplay is the texture to be shown
%drect is the rectangle in which todisplay is shown
% textDown is the distance between the bottom of drect and the text that
% will be drawn on screen


%
%  SP  7/30/08

% commented and adjusted code NW 7/31/08
% switched calls from KeyPress to GetChar
% added tests for ascii codes from useful but non-alpha keypresses
% 8/1/08 NW


    center = [rect(3)/2 rect(4)/2];  %find center coordinates


    Screen('DrawTexture',w,todisplay,[],drect); %draws the  imageabout which the subject is being asked to respond, which is passed in as a parameter of this funciton
    message = 'What is this image? \n Type your response and press enter when you are finished.';%text we would like to show - I have them pressing the space bar instead of 'any key' because pressing enter automatically advances to the next picture
    width = RectWidth(Screen('TextBounds',w,message)); %TextBounds finds the width of that text
    textColor = 0;
    Screen('TextFont', w, 'Arial');
    Screen('TextSize', w, 30);
    DrawFormattedText(w, message, 'center', textDown-.07, 0); %now we can draw it at the center of the screen
    Screen('Flip',w);  %show our text

    FlushEvents('keydown'); %remove keypresses from queue
    
    %     initialize response string
    response = [];

    %main loop to collect responses
    while 1  %run this loop until we escape it using the break command


        lastresponse = response; %keep track of anychanges to response
        
        letter2add = GetChar(0,1);  %need to return ascii code to distinguish return, backspace, etc
        
        if letter2add == 10  %check for return
            break;
        end
        if letter2add == 8%removes the last entered key if the delete key is pressed
            if ~isempty(response) %make sure the array has char to delete
                response = response(1:end-1); %remove last char
                letter2add = '';%nothing will be added to the response array
            end
        end

        response = [response letter2add]; %get the keypress add it to the response

        if isempty(response) %need filler to compute size of text bounds when string is empty
            response = ' ';
        end

        if ~strcmp(response,lastresponse) %only redraw screen if something changed
            width = RectWidth(Screen('TextBounds',w,response)); %TextBounds finds he width of the text in response
            Screen('DrawText',w,response,center(1)-width/2, textDown,0); %now we can draw it at the center of the screen
            Screen('DrawTexture',w,todisplay, [], drect); %draws the imageabout which the subject is being asked to respond, which is passed in as a parameter of this funciton
            Screen('Flip',w);  %show our text
        end
    end


end