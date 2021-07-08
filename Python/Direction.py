import time, serial
from psychopy import visual, core, event
import pandas as pd
import numpy as np
import time

def text(screen):

    text_R3 = visual.TextStim(screen, text=">>>", height=100, color=[1, 1, 1])
    text_R2 = visual.TextStim(screen, text=">>", height=100, color=[1, 1, 1])
    text_R1 = visual.TextStim(screen, text=">", height=100, color=[1, 1, 1])

    text_L3 = visual.TextStim(screen, text="<<<", height=100, color=[1, 1, 1])
    text_L2 = visual.TextStim(screen, text="<<", height=100, color=[1, 1, 1])
    text_L1 = visual.TextStim(screen, text="<", height=100, color=[1, 1, 1])

    text_c = visual.TextStim(screen, text=" + ", height=100, color=[1, 1, 1])

    return text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c


def Direction(tr, original, opposite, screen, port):

    global text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c

    text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c = text(screen)

    # ----------------------------- Psychopy Window & Serial Write ------------------------------#
    origin = globals()['text_{}3'.format(original)]
    opposite = globals()['text_{}3'.format(opposite)]

    if ( tr+1 >= 1 and tr <= 7 ) or (tr+1 >= 15 and tr+1 <= 17) or (tr+1 == 22 or tr+1 ==23 or tr+1 == 26) or (tr+1 == 27 or tr+1 == 30):      # Represent Origin direction
        # Train 1 set               //     Test session 1 - 3       //    Test session 2 - 3                 //  Test sesstion 3 - standard / opposite

        text_c.draw()       # Draw Fixation
        screen.flip()

        # Interval
        time.sleep(0.5)

        origin.draw()      # Draw Direction
        screen.flip()


    elif ( tr+1 >= 8 and tr+1 <= 14 ) or (tr+1 >= 18 and tr+1 <= 20) or (tr+1 == 21 or tr+1 == 24 or tr+1 == 25) or (tr+1 == 28 or tr+1 == 29):      # Represent Opposite direction
        # Train 2 set                //     Test 1 session - 2       //     Test 2 session - 3                  //   Test session 3 - opposite / standard

        text_c.draw()       # Draw Fixation
        screen.flip()

        # Interval
        time.sleep(0.5)

        opposite.draw()      # Draw Direction
        screen.flip()


def switching(tr, check, original, opposi, screen):

    global text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c

    text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c = text(screen)

    origin = globals()['text_{}3'.format(original)]
    opposite = globals()['text_{}3'.format(opposi)]

    if original == 'L':

        if tr+1 == 27:  # 31s   # origin / opposite

            if check == 29:
                origin = globals()['text_{}2'.format(original)]
                origin.draw()
                screen.flip()
            if check == 30:
                origin = globals()['text_{}1'.format(original)]
                origin.draw()
                screen.flip()
            if check == 31:
                opposite.draw()
                screen.flip()

        if tr+1 == 28:  # 32s   # opposite / origin

            if check == 30:
                opposite = globals()['text_{}2'.format(opposi)]
                opposite.draw()
                screen.flip()
            if check == 31:
                opposite = globals()['text_{}1'.format(opposi)]
                opposite.draw()
                screen.flip()
            if check == 32:
                origin.draw()
                screen.flip()

        if tr+1 == 29:  # 34s   # opposite / origin

            if check == 32:
                opposite = globals()['text_{}2'.format(opposi)]
                opposite.draw()
                screen.flip()
            if check == 33:
                opposite = globals()['text_{}1'.format(opposi)]
                opposite.draw()
                screen.flip()
            if check == 34:
                origin.draw()
                screen.flip()

        if tr+1 == 30:  # 28s   # origin / opposite

            if check == 26:
                origin = globals()['text_{}2'.format(original)]
                origin.draw()
                screen.flip()
            if check == 27:
                origin = globals()['text_{}1'.format(original)]
                origin.draw()
                screen.flip()
            if check == 28:
                opposite.draw()
                screen.flip()