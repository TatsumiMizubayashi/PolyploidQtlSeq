﻿using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.Share
{
    /// <summary>
    /// 使用するスレッド数
    /// </summary>
    public class ThreadNumber
    {
        /// <summary>
        /// 使用するスレッド数の最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// 使用するスレッド数の最大値
        /// </summary>
        public const int MAXIMUM = 50;

        /// <summary>
        /// 使用するスレッド数の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 20;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "t";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "thread";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Number of threads to use.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -t option must be an integer greater than or equal to 1 and less than or equal to 50.";

        /// <summary>
        /// 使用するスレッド数を作成する。
        /// </summary>
        /// <param name="number">スレッド数</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public ThreadNumber(int number, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, number, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 使用するスレッド数を作成する。
        /// </summary>
        /// <param name="number">スレッド数</param>
        public ThreadNumber(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }

        /// <summary>
        /// 使用するスレッド数を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns></returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
