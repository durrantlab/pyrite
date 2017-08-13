# DurBlend is a library for easily creating blender plugins.
# Copyright (C) 2017  Jacob D. Durrant
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import bpy
import time

class MessagesClass():
    """
    A class for displaying messages in UI panels. Good for error messages, for
    example.
    """

    messages = {}

    def __init__(self):
        """
        Overwritten by children.
        """
        
        pass

    def send_message(self, id, msg, duration=5.0, msg_otherwise=""):
        """
        Save a message to be displayed in the panel.

        :param str id: An id uniquely identifying this message.

        :param str msg: The message.

        :param float duration: How many seconds the message should be visible.

        :param str msg_otherwise: The message to display if this message isn't
                   being shown. (Not currently implemented, I think...)
        """

        self.messages[id] = {
            "msg": msg,
            "start_time": int(time.time()),
            "duration": duration,
            "msg_otherwise": msg_otherwise  # Not implemented...
        }

    def display_message(self, id, panel):  
        """
        Display a message.

        :param str id: An id uniquely identifying this message.

        :param PanelParentClass panel: The panel object where the message is
                                displayed.
        """

        # Note that must always be called in the context of a panel.
        if id in self.messages:
            msg = self.messages[id]
            if (int(time.time()) - msg["start_time"] < msg["duration"]):
                panel.ui.label(msg["msg"])
            else:
                del self.messages[id]