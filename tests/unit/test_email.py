from itertools import groupby

from mockito import captor, mock, verify, when
from mockito.matchers import ArgumentCaptor
from python_http_client.client import Response
from sendgrid import Mail, SendGridAPIClient

from casp.services.utils.email_utils import EmailSender

BASE_TEMPLATE_DIR: str = "tests/unit/test_email_templates"


class TestEmail:
    """
    Test the email sending and templating functionality.
    """

    def setup_method(self) -> None:
        # Mock the SendGrid client
        self.mock_sendgrid_client = mock(SendGridAPIClient)
        when(self.mock_sendgrid_client).send(...).thenReturn(mock(Response))

    def test_send_email(self) -> None:
        from_email: str = "foo@foo.com"
        email_sender = EmailSender(
            sendgrid_client=self.mock_sendgrid_client, from_address=from_email, template_dir=BASE_TEMPLATE_DIR
        )

        to_email: str = "bar@bar.com"
        subject: str = "Bleep Bloop"
        email_sender._send_email(
            to=to_email, subject=subject, content_path_html="email.html", content_path_plain="email.txt", foo="stuff"
        )

        mail_captor: ArgumentCaptor[Mail] = captor()
        verify(self.mock_sendgrid_client).send(mail_captor)
        mail: Mail = mail_captor.value
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
        email_sender = EmailSender(sendgrid_client=None, from_address=from_email)

        # Neither of these should send an email
        email_sender.send_welcome_email(email="bar@bar.com", key="superdupersecret")
        email_sender.send_new_key_email(email="bar@bar.com", key="superdupersecret")

        verify(self.mock_sendgrid_client, times=0).send(...)
