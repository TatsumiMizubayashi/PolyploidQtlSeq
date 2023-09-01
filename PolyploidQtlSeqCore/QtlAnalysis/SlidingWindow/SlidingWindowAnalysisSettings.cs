namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// スライディングウインドウ解析設定
    /// </summary>
    internal class SlidingWindowAnalysisSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [WindowSize.SHORT_NAME] = WindowSize.LONG_NAME,
            [WindowSize.LONG_NAME] = WindowSize.LONG_NAME,

            [StepSize.SHORT_NAME] = StepSize.LONG_NAME,
            [StepSize.LONG_NAME] = StepSize.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// スライディングウインドウ解析設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public SlidingWindowAnalysisSettings(ISlidingWindowAnalysisSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            WindowSize = new WindowSize(settingValue.WindowSize, parameterDictionary, userOptionDictionary);
            StepSize = new StepSize(settingValue.StepSize, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// Window Size
        /// </summary>
        public WindowSize WindowSize { get; }

        /// <summary>
        /// StepSize
        /// </summary>
        public StepSize StepSize { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                WindowSize.ToParameterFileLine(),
                StepSize.ToParameterFileLine()
            };
        }
    }
}
