using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeq.Options.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// コマンドの基底クラス
    /// </summary>
    [HelpOption("-h|--help", Description = "Show help message.")]
    internal abstract class CommandBase
    {
        /// <summary>
        /// コマンドを実行する。
        /// </summary>
        /// <param name="app">App</param>
        /// <returns>終了コード</returns>
        public abstract Task<int> OnExecuteAsync(CommandLineApplication app);


        /// <summary>
        /// パラメーターファイルから設定を読み込む。
        /// </summary>
        /// <param name="filePath">パラメーターファイルPath</param>
        /// <param name="title">タイトル</param>
        /// <param name="options">オプション</param>
        /// <param name="app">app</param>
        protected static void LoadParameterFile(string filePath, string title, OptionCollection options, CommandLineApplication app)
        {
            if (string.IsNullOrEmpty(filePath)) return;

            Console.WriteLine($"Load settings from {filePath}.");
            var parameterFile = new ParameterFile(title, options);
            var paramsDict = parameterFile.Parse(filePath);
            options.SetValues(paramsDict, app.Options);
        }

        /// <summary>
        /// データ検証を行う。
        /// </summary>
        /// <param name="options">オプション</param>
        /// <returns>エラーがある場合はtrue</returns>
        protected static bool Validation(OptionCollection options)
        {
            var errorValidations = options.Validation();
            if (errorValidations.Length == 0) return false;

            foreach (var error in errorValidations)
            {
                error.Print();
            }

            return true;
        }

        /// <summary>
        /// パラメーターファイルを作成する。
        /// </summary>
        /// <param name="filePath">パラメーターファイルPath</param>
        /// <param name="title">タイトル</param>
        /// <param name="options">オプション</param>
        protected static void CreateParameterFile(string filePath, string title, OptionCollection options)
        {
            var parameterFile = new ParameterFile(title, options);
            parameterFile.Create(filePath);
        }
    }
}
