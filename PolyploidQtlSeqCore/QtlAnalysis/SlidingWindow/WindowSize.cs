using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Windowサイズ
    /// </summary>
    public class WindowSize
    {
        /// <summary>
        /// window sizeの最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// window sizeの最大値
        /// </summary>
        public const int MAXIMUM = 100000;

        /// <summary>
        /// window sizeの規定値
        /// </summary>
        public const int DEFAULT = 100;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "w";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "window";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Window size (kbp) of the sliding window analysis.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -w option must be an integer greater than 1 and less than or equal to 100,000.";

        /// <summary>
        /// WindowSizeを作成する。
        /// </summary>
        /// <param name="size">サイズ(kbp)</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public WindowSize(int size, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            KbpValue = OptionValue.GetValue(LONG_NAME, size, parameterDictionary, userOptionDictionary);

            if (KbpValue < MINIMUM || KbpValue > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);

            BpValue = KbpValue * 1000;
        }

        /// <summary>
        /// WindowSize(kbp)を取得する。
        /// </summary>
        public int KbpValue { get; }

        /// <summary>
        /// WindowSize(bp)を取得する。
        /// </summary>
        public int BpValue { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{KbpValue}";
        }
    }
}
