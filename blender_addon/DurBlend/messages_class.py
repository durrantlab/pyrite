import bpy
import time

class MessagesClass():
    messages = {}

    def __init__(self):
        # A place to store messages that can be passed
        #self.messages = {}
        pass

    def send_message(self, id, msg, duration=5.0, msg_otherwise=""):
        self.messages[id] = {
            "msg": msg,
            "start_time": int(time.time()),
            "duration": duration,
            "msg_otherwise": msg_otherwise  # Not implemented...
        }

    def display_message(self, id, panel):  # Note that must always be called in the context of a panel.
        if id in self.messages:
            msg = self.messages[id]
            if (int(time.time()) - msg["start_time"] < msg["duration"]):
                panel.ui.label(msg["msg"])
            else:
                del self.messages[id]