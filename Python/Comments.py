from psychopy import visual, core, event
import pandas as pd
import time

correct = []
answer = []

def Question(j, file, path, screen):

    file = pd.read_excel(path + "/AAD/Python/question.xlsx")

    try :
        # Question 1
        text3 = visual.TextStim(screen, text = file.tweenty_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 2
        text3 = visual.TextStim(screen, text = file.tweenty_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A2[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 3
        text3 = visual.TextStim(screen, text = file.journey_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.journey_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 4
        text3 = visual.TextStim(screen, text = file.journey_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
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
    time.sleep(3)

    return correct, answer



def Comments(tr,path, screen):
    file_2 = pd.read_excel(path + "/AAD/Python/Comments.xlsx")

    while True:
        try:
            # Comment
            for i in range(0, 11):

                if tr == 'intro':
                    print("Intro & Train Session 1")
                    text = visual.TextStim(screen, text=file_2.intro[i], height=50, color=[1, 1, 1], wrapWidth=2000)

                elif tr+1 == 7:     # Train Session 2
                    print("Train Seession 2")
                    text = visual.TextStim(screen, text=file_2.train2[i], height=50, color=[1, 1, 1], wrapWidth=2000)

                elif tr+1 == 14:    # Test session 1
                    print("Test session 1")
                    text = visual.TextStim(screen, text=file_2.test1[i], height=50, color=[1, 1, 1], wrapWidth=2000)

                elif tr+1 == 20:    # Test session 2
                    print("Test session 2")
                    text = visual.TextStim(screen, text=file_2.test2[i], height=50, color=[1, 1, 1], wrapWidth=2000)

                elif tr+1 == 26:    # Test session 3
                    print("Test session 3")
                    text = visual.TextStim(screen, text=file_2.test3[i], height=50, color=[1, 1, 1], wrapWidth=2000)

                text.draw()
                screen.flip()

                key = event.waitKeys(keyList=["space", "escape"], clearEvents=True)
                if key == ["escape"]:
                    core.quit()

        except:
            pass
            break