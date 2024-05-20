from itertools import groupby

import sendgrid
from mockito import captor, mock, unstub, verify, when
from mockito.matchers import ArgumentCaptor
from python_http_client.client import Response

from casp.services.utils.email_utils import EmailSender

BASE_TEMPLATE_DIR: str = "tests/unit/test_email_templates"


class TestEmail:
    """
    Test the email sending and templating functionality.
    """

    def setup_method(self) -> None:
        # Mock the SendGrid client
        self.mock_sendgrid_client = mock(sendgrid.SendGridAPIClient)
        when(self.mock_sendgrid_client).send(...).thenReturn(mock(Response))
        when(sendgrid).SendGridAPIClient(...).thenReturn(self.mock_sendgrid_client)

    def teardown_method(self) -> None:
        unstub()

    def test_send_email(self) -> None:
        from_email: str = "foo@foo.com"
        sendgrid_key: str = "super_secret_key"
        email_sender = EmailSender(sendgrid_key=sendgrid_key, from_address=from_email, template_dir=BASE_TEMPLATE_DIR)

        verify(sendgrid, times=1).SendGridAPIClient(sendgrid_key)

        to_email: str = "bar@bar.com"
        subject: str = "Bleep Bloop"
        email_sender._send_email(
            to=to_email,
            subject=subject,
            content_path_html="email.html",
            content_path_plain="email.txt",
            sub_values={"foo": "stuff"},
        )

        mail_captor: ArgumentCaptor[sendgrid.Mail] = captor()
        verify(self.mock_sendgrid_client).send(mail_captor)
        mail: sendgrid.Mail = mail_captor.value
        assert mail.from_email.email == from_email
        assert mail.personalizations[0].tos[0]["email"] == to_email
        assert mail.subject.subject == subject
        # Verify that the email content is correct
        assert len(mail.contents) == 2
        for mime_type, content in groupby(mail.contents, lambda c: c.mime_type):
            content = list(content)
            if mime_type == "text/plain":
                assert len(content) == 1
                assert content[0].content == "Foo is: stuff"
            elif mime_type == "text/html":
                assert len(content) == 1
                assert content[0].content == "Foo is: <div>stuff</div>"

    def test_on_not_send_on_no_client(self) -> None:
        from_email: str = "foo@foo.com"
        email_sender = EmailSender(sendgrid_key=None, from_address=from_email)

        verify(sendgrid, times=0).SendGridAPIClient(...)

        # Neither of these should send an email
        email_sender.send_welcome_email(email="bar@bar.com", key="superdupersecret")
        email_sender.send_new_key_email(email="bar@bar.com", key="superdupersecret")

        verify(self.mock_sendgrid_client, times=0).send(...)
