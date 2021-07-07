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


def Direction(tr, w, original, opposite, screen, port):

    global text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c

    text_R3, text_R2, text_R1, text_L3, text_L2, text_L1, text_c = text(screen)

    # ----------------------------- Psychopy Window & Serial Write ------------------------------#
    origin = globals()['text_{}3'.format(original)]
    opposite = globals()['text_{}3'.format(opposite)]

    # Press Button for start
    if tr+1 == 1:

        # Draw Fixation
        text_c.draw()
        screen.flip()

        key2 = event.getKeys()
        if key2 == ["space"]:
            # Ready
            time.sleep(3)
            port.write(b'1')

            origin.draw()
            screen.flip()
            w = 0

    elif ( tr+1 >= 2 and tr <= 7 ) or (tr+1 >= 15 and tr+1 <= 17) or (tr+1 == 22 or tr+1 ==23 or tr+1 == 26) or (tr+1 == 27 or tr+1 == 30):      # Represent Origin direction
        # Train 1 set               //     Test session 1 - 3       //    Test session 2 - 3                 //  Test sesstion 3 - standard / opposite

        # Interval
        time.sleep(0.5)

        text_c.draw()       # Draw Fixation
        screen.flip()
        port.write(b'1')    # To arduino

        origin.draw()      # Draw Direction
        screen.flip()
        w = 0

    elif ( tr+1 >= 8 and tr+1 <= 14 ) or (tr+1 >= 18 and tr+1 <= 20) or (tr+1 == 21 or tr+1 == 24 or tr+1 == 25) or (tr+1 == 28 or tr+1 == 29):      # Represent Opposite direction
        # Train 2 set                //     Test 1 session - 2       //     Test 2 session - 3                  //   Test session 3 - opposite / standard

        # Interval
        time.sleep(0.5)

        text_c.draw()       # Draw Fixation
        screen.flip()
        port.write(b'1')    # To arduino

        opposite.draw()      # Draw Direction
        screen.flip()
        w = 0

    return w


def switching(tr, i, original, opposite, screen):

    text_R3 = visual.TextStim(screen, text=">>>", height=100, color=[1, 1, 1])
    text_R2 = visual.TextStim(screen, text=">>", height=100, color=[1, 1, 1])
    text_R1 = visual.TextStim(screen, text=">", height=100, color=[1, 1, 1])

    text_L3 = visual.TextStim(screen, text="<<<", height=100, color=[1, 1, 1])
    text_L2 = visual.TextStim(screen, text="<<", height=100, color=[1, 1, 1])
    text_L1 = visual.TextStim(screen, text="<", height=100, color=[1, 1, 1])

    text_c = visual.TextStim(screen, text=" + ", height=100, color=[1, 1, 1])

    origin = globals()['text_{}3'.format(original)]
    opposite = globals()['text_{}3'.format(opposite)]

    if original == 'L':

        if tr == 27:  # 31s   # origin / opposite

            if i == 29:
                origin = globals()['text_{}2'.format(original)]
                origin.draw()
                screen.flip()
            if i == 30:
                origin = globals()['text_{}1'.format(original)]
                origin.draw()
                screen.flip()
            if i == 31:
                opposite.draw()
                screen.flip()

        if tr == 28:  # 29s   # opposite / origin

            if i == 27:
                opposite = globals()['text_{}2'.format(opposite)]
                opposite.draw()
                screen.flip()
            if i == 28:
                opposite = globals()['text_{}1'.format(opposite)]
                opposite.draw()
                screen.flip()
            if i == 29:
                origin.draw()
                screen.flip()

        if tr == 29:  # 30s   # opposite / origin

            if i == 28:
                opposite = globals()['text_{}2'.format(opposite)]
                opposite.draw()
                screen.flip()
            if i == 29:
                opposite = globals()['text_{}1'.format(opposite)]
                opposite.draw()
                screen.flip()
            if i == 30:
                origin.draw()
                screen.flip()

        if tr == 30:  # 28s   # origin / opposite

            if i == 26:
                origin = globals()['text_{}2'.format(original)]
                origin.draw()
                screen.flip()
            if i == 27:
                origin = globals()['text_{}1'.format(original)]
                origin.draw()
                screen.flip()
            if i == 28:
                opposite.draw()
                screen.flip()