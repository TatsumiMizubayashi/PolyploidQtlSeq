using Kurukuru;
using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.QualityControl;

namespace PolyploidQtlSeqCore.Application.QualityControl
{
    /// <summary>
    /// FastpによるQuality Control
    /// </summary>
    public class FastpQualityControl
    {
        private const string LOG_DIR_NAME = "Log";

        //private readonly FastpQualityControlCommandOptions _qcCommandOption;
        private readonly FastpQualityControlSettings _settings;

        /// <summary>
        /// Fastp Quality Controlを作成する。
        /// </summary>
        /// <param name="optionValues">コマンドオプションの値</param>
        /// <param name="options">オプションリスト</param>
        public FastpQualityControl(IFastpQualityControlCommandOptions optionValues,
            IReadOnlyCollection<CommandOption> options)
        {
            //_qcCommandOption = new FastpQualityControlCommandOptions(optionValues, options);
        }

        /// <summary>
        /// Fastp Quality Controlインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public FastpQualityControl(IFastpQualityControlSettingValue settingValue)
        {
            _settings = new FastpQualityControlSettings(settingValue);
        }

        /// <summary>
        /// Quality Controlを実行する。
        /// </summary>
        /// <returns>終了コード</returns>
        public async ValueTask<int> RunAsync()
        {
            var inputFastqFilePairs = _settings.InputRawFastqDirectory.ToFastqFilePairs();
            var fastpCommonOption = _settings.ToFastpCommonOption();

            var outputDir = fastpCommonOption.OutputDirectory;
            var logDirPath = outputDir.CreateSubDir(LOG_DIR_NAME);

            var code = 0;
            Log.Clear();
            CommandLog.Clear();

            foreach (var inputFilePair in inputFastqFilePairs)
            {
                await Spinner.StartAsync($"{inputFilePair.BaseName} quality control ...", async spinner =>
                {
                    try
                    {
                        await Fastp.RunAsync(inputFilePair, fastpCommonOption);
                        spinner.Succeed($"{inputFilePair.BaseName} completed");
                    }
                    catch(Exception ex)
                    {
                        // 例外が出ても次のサンプルを処理する。
                        code = 1;
                        spinner.Fail($"{inputFilePair.BaseName} error");
                        Console.Error.WriteLine(ex.Message);
                        Console.Error.WriteLine(ex.StackTrace);
                    }
                    finally
                    {
                        var logFilePath = Path.Combine(logDirPath, inputFilePair.BaseName + ".Log.txt");
                        Log.Save(logFilePath);
                        Log.Clear();
                    }
                });

            }

            var commandListFilePath = outputDir.CreateFilePath(LOG_DIR_NAME, "QC Command List.txt");
            CommandLog.Save(commandListFilePath);

            return code;
        }
    }
}
