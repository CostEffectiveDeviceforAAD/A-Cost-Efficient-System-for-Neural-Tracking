from psychopy import visual, core, event
import pandas as pd
import time
import numpy as np

correct = []
answer = []

def practice(p, path, screen):

    print("practice")
    file = pd.read_excel(path + "/AAD/Python/prePractice.xlsx")

    text = visual.TextStim(screen, text = "<<<", height=100, color=[1, 1, 1], wrapWidth=1500)
    text.draw()
    screen.flip()
    time.sleep(15.3)

    # Question 1
    text = visual.TextStim(screen, text = file.tweenty_Q1_practice[p], height=50, color=[1, 1, 1], wrapWidth=1500)
    text.draw()
    screen.flip()

    event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    # Question 2
    text = visual.TextStim(screen, text = file.journey_Q1_practice[p], height=50, color=[1, 1, 1], wrapWidth=1500)
    text.draw()
    screen.flip()

    event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    text = visual.TextStim(screen, text = "+", height=100, color=[1, 1, 1], wrapWidth=1500)
    text.draw()
    screen.flip()
    time.sleep(0.5)


def Question(j, path, screen):
    correct = []
    file = pd.read_excel(path + "/AAD/Python/question.xlsx")

    try :
        # Question 1
        text3 = visual.TextStim(screen, text = file.tweenty_Q1[j], height=50, color=[1, 1, 1], wrapWidth=1500)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 2
        text3 = visual.TextStim(screen, text = file.tweenty_Q2[j], height=50, color=[1, 1, 1], wrapWidth=1500)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A2[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 3
        text3 = visual.TextStim(screen, text = file.journey_Q1[j], height=50, color=[1, 1, 1], wrapWidth=1500)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.journey_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 4
        text3 = visual.TextStim(screen, text = file.journey_Q2[j], height=50, color=[1, 1, 1], wrapWidth=1500)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.journey_A2[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

    except:
        correct.append("N")

    # Interval
    text = visual.TextStim(screen, text=" + ", height=100, color=[1, 1, 1], wrapWidth=2000)
    text.draw()
    screen.flip()
    #time.sleep(3)

    return correct, answer



def Comments(tr,path, screen):
    n = []
    i=0
    file_2 = pd.read_excel(path + "/AAD/Python/Comments.xlsx")

    try:
        # Comment
        for i in range(0,15):

            if tr == 'intro':
                print("Intro")
                text = visual.TextStim(screen, text=file_2.intro[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.intro[i]

            elif tr+1 == 1:     # Train Session 1
                print("Train Seession 1  ")
                text = visual.TextStim(screen, text=file_2.train1[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.train1[i]

            elif tr == 7:     # Train Session 2
                print("Train Seession 2")
                text = visual.TextStim(screen, text=file_2.train2[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.train2[i]

            elif tr == 14:    # Test session 1
                print("Test session 1")
                text = visual.TextStim(screen, text=file_2.test1[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.test1[i]

            elif tr == 20:    # Test session 2
                print("Test session 2")
                text = visual.TextStim(screen, text=file_2.test2[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.test2[i]

            elif tr == 26:    # Test session 3
                print("Test session 3")
                text = visual.TextStim(screen, text=file_2.test3[i], height=37, color=[1, 1, 1], wrapWidth=1500)
                n = file_2.test3[i]

            try:
                np.isnan(n)
                break
            except:
                text.draw()
                screen.flip()

            key = event.waitKeys(keyList=["space", "escape"], clearEvents=True)
            if key == ["escape"]:
                core.quit()
    except:
        pass

