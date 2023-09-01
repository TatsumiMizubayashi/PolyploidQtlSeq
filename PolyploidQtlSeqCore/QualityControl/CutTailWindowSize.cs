using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 3'末端トリム時のウインドウサイズ
    /// </summary>
    public class CutTailWindowSize
    {
        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 1;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "W";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "cutTailWindowSize";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Window size when trimmed at 3' end. cut_tail_window_size in fastp.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -W option must be an integer greater than or equal to 1 and less than or equal to 100.";

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを作成する。
        /// </summary>
        /// <param name="windowSize">ウインドウサイズ</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public CutTailWindowSize(int windowSize, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, windowSize, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを作成する。
        /// </summary>
        /// <param name="windowSize">ウインドウサイズ</param>
        public CutTailWindowSize(int windowSize)
        {
            if (windowSize < MINIMUM || windowSize > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(windowSize));

            Value = windowSize;
        }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns></returns>
        internal string ToFastpArg()
        {
            return $"--cut_tail_window_size {Value}";
        }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }
    }
}
