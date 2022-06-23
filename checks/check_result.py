import enum
# Order important - later statuses have higher precedence.
CheckStatus = enum.Enum('CheckStatus', ['GREEN', 'YELLOW', 'RED'])


class CheckResult:
    def __init__(self, status, message, batch, context, table_type):
        self.status = status
        self.context = context or 'General'
        self.table_type = table_type or 'General'
        if batch is None:
            self.message = message
        else:
            self.message = "["+batch+"] "+str(message)
    def __repr__(self):
        return '%s: %s' % (self.status, self.message)
    def to_json(self):
        return {
            'heading': self.table_type,
            'tabletype': self.context,
            'color': self.status,
            'message': self.message
            }

def green(message, batch, context, table_type):
    return CheckResult('Green', message, batch, context, table_type)


def yellow(message, batch, context, table_type):
    return CheckResult('Yellow', message, batch, context, table_type)


def red(message, batch, context, table_type):
    return CheckResult('Red', message, batch, context, table_type)
